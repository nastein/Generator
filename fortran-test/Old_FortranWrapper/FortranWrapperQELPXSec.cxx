//____________________________________________________________________________
/*
 Copyright (c) 2003-2021, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Noah Steinberg <nsteinbe \at fnal.gov>
         Syrian Truong <struong@fnal.gov>
         Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Acclerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#include <algorithm>

#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearUtils.h"

#include "FortranWrapperQELPXSec.h"
#include "FortranWrapperXSecIntegrator.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
extern"C"
{
void diracmatrices_(double *xmn_in);
}

//____________________________________________________________________________
extern"C"
{
void compute_hadron_tensor_(double *xq, double *w, double *wt, double *xk, double *xp, std::complex<double> resp[4][4]/*response tensor*/, double *nuphi, double *f1v, double *f2v, double *ffa);
}

//____________________________________________________________________________
FortranWrapperQELPXSec::FortranWrapperQELPXSec() :
XSecAlgorithmI("genie::FortranWrapperQELPXSec")
{

}
//____________________________________________________________________________
FortranWrapperQELPXSec::FortranWrapperQELPXSec(string config) :
XSecAlgorithmI("genie::FortranWrapperQELPXSec", config)
{

}
//____________________________________________________________________________
FortranWrapperQELPXSec::~FortranWrapperQELPXSec()
{

}
//____________________________________________________________________________
double FortranWrapperQELPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{

  if ( !this->ValidProcess(interaction)) return 0.;
  if ( !this->ValidKinematics(interaction) ) return 0.;

  // Get kinematics and init-state parameters
  // Particles are stored in the lab frame
  const Kinematics&   kinematics = interaction->Kine();
  const InitialState& init_state = interaction->InitState();
  const Target& target = init_state.Tgt();

  //Final state lepton
  TLorentzVector fsLeptonP4 = kinematics.FSLeptonP4();
  TVector3 leptonMom3 = fsLeptonP4.Vect();
  double E_lep = fsLeptonP4.E();
  double P_lep = fsLeptonP4.P();

  //Final state nucleon 
  TLorentzVector p4Nf = kinematics.HadSystP4();
  TVector3 NfMom3 = p4Nf.Vect();
  double E_Nf = p4Nf.E();
  double P_Nf = p4Nf.P();

  //Initial state lepton
  TLorentzVector* tempProbeP4 = init_state.GetProbeP4( kRfLab );
  TLorentzVector probeP4 = *tempProbeP4;
  delete tempProbeP4;
  TVector3 probeMom3 = probeP4.Vect();
  double p_probe = probeP4.P();

  //Initial state nucleon
  TLorentzVector p4Ni = target.HitNucP4();
  double mNi = target.HitNucMass(); 
  double pNi = p4Ni.P();
  double ENi = p4Ni.E();
  double E_NiOnShell = std::sqrt(pNi*pNi + mNi*mNi);
  //Make On shell 4 vector for nucleon
  TLorentzVector p4NiOnShell = TLorentzVector(p4Ni.Vect(), E_NiOnShell);
  TVector3 NiMom3 = p4Ni.Vect(); 
  double epsilon_B = E_NiOnShell - p4Ni.E();

  // Check outgoing nucleon momentum against kF. If it is Pauli-blocked, then
  // just return zero without bothering to do the rest of the calculation.
  double pNf = p4Nf.P();
  double kF = fNuclModel->LocalFermiMomentum( target, target.HitNucPdg(), pNi );
  if ( pNf < kF ) return 0.;

  // Noemi's fortran wrapper requires that q3Vec is in the Z direction, 
  // so rotate the lepton momenta
  TVector3 q3Vec = probeMom3 - leptonMom3;
  TVector3 zvec(0.0, 0.0, 1.0);
  TVector3 rot = (q3Vec.Cross(zvec) ).Unit(); //Vector to rotate about
  // Angle between the z direction and q
  double angle = zvec.Angle( q3Vec );

  // Handle the edge case where q3Vec is along -z so the cross product above vanishes
  if ( q3Vec.Perp() == 0. && q3Vec.Z() < 0. ) {
    rot = TVector3( 0., 1., 0.);
    angle = kPi;
  }

  if ( rot.Mag() >= kASmallNum ) {
    probeMom3.Rotate(angle, rot);
    probeP4.SetVect(probeMom3);
    
    leptonMom3.Rotate(angle, rot);
    fsLeptonP4.SetVect(leptonMom3);

    NiMom3.Rotate(angle,rot);
    p4Ni.SetVect(NiMom3);
    p4NiOnShell.SetVect(NiMom3);

    NfMom3.Rotate(angle,rot);
    p4Nf.SetVect(NfMom3);     
  }

  // 4-momentum transfer, qTildeP4 corresponds to the de Forest prescription for handling off-shell initial struck neclon
  TLorentzVector qP4 = probeP4 - fsLeptonP4;
  TLorentzVector qTildeP4 = p4Nf - p4NiOnShell;

  double Q2 = -1. * qP4.M2();
  double Q2tilde = -1. * qTildeP4.M2();

  // Apply the global Q^2 cut imposed by GENIE
  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if ( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics::electromagnetic::kMinQ2Limit; // EM limit
  //std::cout << "Q2min = " << Q2min << "\n";
  //Q2min = 0.164;

  if ( Q2 < Q2min ) return 0.;

  //If the binding energy correction causes an unphysical value of q0Tilde or Q2Tilde, just return 0.
  if ( qTildeP4.E() <= 0. && init_state.Tgt().IsNucleus() &&
       !interaction->TestBit(kIAssumeFreeNucleon) ) return 0.;
  if ( Q2tilde <= 0. ) return 0.;

  // Cross section 
  //double xsec = 1.0;
  double xsec = 1./(64. * kPi * kPi * E_lep * E_Nf * p_probe * ENi);

  

  // Scale cross section by the number of active nucleons
  int hit_nuc_pdg = target.HitNucPdg();
  //Number of active nucleons
  int num_active = pdg::IsProton(hit_nuc_pdg) ? target.Z() : target.N(); 

  xsec *= num_active;

  //Get the correct coupling and form factors model for the current interaction
  double coupling_factor = 1.;
  const genie::ProcessInfo& proc_info = interaction->ProcInfo();

  if ( proc_info.IsEM() ) {
    fFormFactors.SetModel( fEMFormFactorsModel );
   // coupling_factor = kAem2 * E_lep/ (probeP4.E() * Q2*Q2);
  }

  else if ( proc_info.IsWeakNC() ) {
    fFormFactors.SetModel( fNCFormFactorsModel );
    coupling_factor = kGF2 / 2.;
  }

  else if ( proc_info.IsWeakCC() ) {
    fFormFactors.SetModel( fCCFormFactorsModel );
    coupling_factor = kGF2 * fCos8c2 * P_lep/ (16.0 * kPi * kPi * p_probe) ;
  }

  xsec *= coupling_factor;
  //Set Q2 to Q2tilde while computing form factors
  interaction->KinePtr()->SetQ2( Q2tilde );

  // Evaluate the form factors
  fFormFactors.Calculate( interaction );

  double f1v, f2v, ffa;

  // Let's get the form factors 
  f1v = fFormFactors.F1V();
  f2v = fFormFactors.xiF2V();
  ffa = fFormFactors.FA();
  // Now that we've calculated them, store the true Q2 value
  interaction->KinePtr()->SetQ2( Q2 );

  // Compute the leptonic tensor
  InteractionType_t type = interaction->ProcInfo().InteractionTypeId();
  // Since we rotated the momentum vectors of the probe and resulting lepton
  // We have to explicitly provide LeptonTensor with the 4-vectors we modified
  // or else it will use the lab frame vectors
 // LeptonTensor L_munu(probeP4, fsLeptonP4, init_state.ProbePdg(), type);

  // xmn_in (GeV)
  double pm_GeV = genie::PDGLibrary::Instance()->Find(2212)->Mass(); 
  double nm_GeV = genie::PDGLibrary::Instance()->Find(2112)->Mass(); 
  double xmn_in = 0.5 * (pm_GeV + nm_GeV);

  //q0 (GeV)
  double w = qP4.E();

  //qtilde0 (GeV)
  double wt = qTildeP4.E();
  
  // 3-momentum transfer (GeV)
  double xq = qP4.P();

  // Initial nucleon momentum (GeV)
  double xk = pNi;

  // Final nucleon momentum (GeV)
  double xp = pNf; 

  // Final lepton angle 
  double theta = std::abs(fsLeptonP4.Vect().Angle(probeP4.Vect()));
  
  // Final nucleon phi 
  double nuphi = p4Nf.Phi();

  //Set up the 4x4 hadronic response tensor
  std::complex<double> HadronTensor[4][4]; 
   
  // Set up dirac matrices via calling diracmatrices Fortran90 subroutine
  diracmatrices_(&xmn_in);

  // Calculate hadronic response tensor using Fortran90 subroutine
  // Function of incoming and outgoing nucleon momenta, 
  // outoging lepton angle, energy transfer,  and form factors
  // Hadron Tensor is 4x4 complex<double> returned in (x,y,z,t) convention
 // compute_hadron_tensor_(&xq, &w, &wt, &xk, &xp, HadronTensor, &nuphi, &f1v, &f2v, &ffa);
    
  // We're going to zero out certain components to isolate the transverse/longitudinal parts
  bool trans = false;
  bool longit = false;
  bool apply_current_conservation = false;
    
  /*
    if (apply_current_conservation == true) {
      // These lines ensure gauge invariance, i.e. that the Hadronic current is conserved.
      HadronTensorCons[2][3] = (w/xq)*HadronTensorCons[3][3];
      HadronTensorCons[3][2] = (w/xq)*HadronTensorCons[3][3];
      HadronTensorCons[2][2] = pow(w/xq,2)*HadronTensorCons[3][3];
    } 
   
    
    if (longit) {
      for(int i = 0; i <4 ; i++) {
      for(int j = 0; j < 4; j++) {
        if((i == 2 || i == 3) && (j == 2 || j == 3))
        {}
        else HadronTensor[i][j] = 0.0;
      }
      }
    }

    if (trans) {
      for(int i = 0; i <4 ; i++) {
      for(int j = 0; j < 4; j++) {
        if((i == 0 ||i == 1) && (j == 0 || j ==1))
        {}
        else {
          //std::cout << "i = " << i << ", j = " << j << "\n";
          HadronTensor[i][j] = 0.0;
        }
      }
      }
    }
  
  // Convert tensor to genie::Rank2LorentzTensor
  ManualResponseTensor W_munu(HadronTensor);
    
  //Contract tensor components
  std::complex<double> contraction = L_munu * W_munu;
  if ( std::abs(contraction.imag()) > kASmallNum ) {
     LOG("FortranWrapperQELPXsec", pWARN) << "Tensor contraction has nonvanishing"
      << " imaginary part!";
  }
    
  xsec *= contraction.real();
     
  // Transform to kPSFullNBody phase space
  xsec *= 4 * E_Nf/E_lep;
*/
   //xsec *= 2. * E_Nf / ( kPi * P_Nf * P_Nf * P_lep );
   return xsec;
 
}
//____________________________________________________________________________
double FortranWrapperQELPXSec::Integral(const Interaction* in) const
{
  // We'll take care of this later. We need to use this function to get
  // total cross sections during spline generation.
  double integ = fXSecIntegrator->Integrate( this, in );
  return integ;
}
//____________________________________________________________________________
bool FortranWrapperQELPXSec::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const InitialState& init_state = interaction->InitState();
  const ProcessInfo&  proc_info  = interaction->ProcInfo();

  if ( !proc_info.IsQuasiElastic() ) return false;

  if ( proc_info.IsEM() || proc_info.IsWeakNC() ) return true;

  // For weak CC interactions, check that the hit nucleon "works"
  else if ( proc_info.IsWeakCC() ) {
    int nucleon_pdg = init_state.Tgt().HitNucPdg();
    int probe_pdg = init_state.ProbePdg();

    if ( pdg::IsNeutron(nucleon_pdg) && pdg::IsNeutrino(probe_pdg) ) {
      return true;
    }
    else if ( pdg::IsProton(nucleon_pdg) && pdg::IsAntiNeutrino(probe_pdg) ) {
      return true;
    }
    else return false;
  }

  return true;
}
//____________________________________________________________________________
void FortranWrapperQELPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperQELPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperQELPXSec::LoadConfig(void)
{

  double thc;
  GetParam( "CabibboAngle", thc);
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

 // fphase_space_only = false;
  // Bool to control integration over only phase space
 // GetParam( "PhaseSpaceOnly", fphase_space_only);

 // fphase_space_only = false;

  // Get access to the nuclear model for use in Pauli blocking, etc.
  fNuclModel = dynamic_cast< const NuclearModelI* >(
    this->SubAlg("NuclearModel") );
  assert( fNuclModel );

  // Load XSec Integrator
  fXSecIntegrator = dynamic_cast< const XSecIntegratorI* >(
    this->SubAlg("XSec-Integrator") );
  assert( fXSecIntegrator );

  // load QEL form factors models
  fEMFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("EMFormFactorsAlg") );
  assert( fEMFormFactorsModel );

  fCCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("CCFormFactorsAlg") );
  assert( fCCFormFactorsModel );

  fNCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("NCFormFactorsAlg") );
  assert( fNCFormFactorsModel );
  
  //Attach CC model for now. This will be updated later
  fFormFactors.SetModel( fCCFormFactorsModel );


}
