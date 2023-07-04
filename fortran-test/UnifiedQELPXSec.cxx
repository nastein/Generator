//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Acclerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

#include "UnifiedQELPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
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

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
extern"C"
{
void diracmatrices_(double *xmn_in);
}

////____________________________________________________________________________
extern"C"
{

void compute_hadron_tensor_(double *wt, double *xk_x, double *xk_y, double *xk_z, double *xp_x, double *xp_y, double *xp_z, double *f1v, double *f2v, double *ffa, double *ffp, std::complex<double> resp[4][4]);

}

//____________________________________________________________________________
UnifiedQELPXSec::UnifiedQELPXSec() :
XSecAlgorithmI("genie::UnifiedQELPXSec")
{

}
//____________________________________________________________________________
UnifiedQELPXSec::UnifiedQELPXSec(string config) :
XSecAlgorithmI("genie::UnifiedQELPXSec", config)
{

}
//____________________________________________________________________________
UnifiedQELPXSec::~UnifiedQELPXSec()
{

}
//____________________________________________________________________________
double UnifiedQELPXSec::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  if ( !this->ValidProcess(interaction) ) return 0.;
  if ( !this->ValidKinematics(interaction) ) return 0.;

  if ( kps == kPSQ2fE ) return this->FreeNucXSec( interaction, kps );

  // Get kinematics and init-state parameters
  const Kinematics&   kinematics = interaction->Kine();
  const InitialState& init_state = interaction->InitState();
  const Target& target = init_state.Tgt();

  TLorentzVector* temp_probeP4 = init_state.GetProbeP4( kRfLab );
  TLorentzVector probeP4 = *temp_probeP4;
  delete temp_probeP4;
  double p_probe = probeP4.P();

  TLorentzVector p4Ni = target.HitNucP4();
  double mNi = target.HitNucMass(); // on-shell initial hit nucleon mass
  double pNi = p4Ni.P();
  double E_NiOnShell = std::sqrt(pNi*pNi + mNi*mNi);
  TLorentzVector p4NiOnShell = TLorentzVector(p4Ni.Vect(), E_NiOnShell);
  double epsilon_B = E_NiOnShell - p4Ni.E();

  TLorentzVector lepP4 = kinematics.FSLeptonP4();
  double E_lep = lepP4.E();
  double P_lep = lepP4.P();

  TLorentzVector p4Nf = kinematics.HadSystP4();
  double E_Nf = p4Nf.E();
  double P_Nf = p4Nf.P();
  
  // in the dOmegadE' phase space, but QELEvgen phase space is
  // defined for a conversion between d3k'(Lab)->dcosthetadphi (COM frame), so we have to
  // convert dOmegadE'= d3k'/E'k'

  double xsec = 1.0;
 
  bool testing_phase_space = false;  
  
  //Phase space volume
  if(testing_phase_space) xsec *= fXSecScale / 16./ kPi / kPi / E_lep / E_Nf;

  if(!testing_phase_space) xsec *= fXSecScale / (E_lep * P_lep);

  // If we're dealing with a nuclear target, then apply Pauli blocking as
  // needed  
  if ( target.IsNucleus() && !interaction->TestBit(kIAssumeFreeNucleon) ) {
    double kF = fPauliBlocker->GetFermiMomentum(target, interaction->RecoilNucleonPdg(), target.HitNucPosition());
    //std::cout << "kF = " << kF << "\n";
    if ( p4Nf.P() < kF ) {return 0.;}
  }

  // Scale cross section by the number of active nucleons
  int hit_nuc_pdg = target.HitNucPdg();
  // Number of active nucleons in the target
  int num_active = pdg::IsProton(hit_nuc_pdg) ? target.Z() : target.N();

  xsec *= num_active;

  // Compute form factors using Q2tilde (the effective Q2 value after
  // binding energy corrections)
  TLorentzVector qP4 = probeP4 - lepP4;
  TLorentzVector qTildeP4 = qP4;
  qTildeP4.SetE( qP4.E() - epsilon_B );

  double Q2 = -1. * qP4.M2();
  double Q2tilde = -1. * qTildeP4.M2();

  double Q2min = genie::controls::kMinQ2Limit; // CC/NC limit
  if( interaction->ProcInfo().IsEM() ) Q2min = genie::utils::kinematics::electromagnetic::kMinQ2Limit; // EM limit
  // Make sure Q2 is physical
  if ( Q2 < Q2min ) {
    return 0.;
  }

  // Get the correct couplings and form factors model for the current
  // interaction
  double coupling_factor = 1.;
  const genie::ProcessInfo& proc_info = interaction->ProcInfo();
  if ( proc_info.IsWeakCC() ) {
    coupling_factor = kGF2 * fCos8c2 * P_lep / (32 * kPi * kPi * p_probe);
    fFormFactors.SetModel( fCCFormFactorsModel );
  }

  else if ( proc_info.IsWeakNC() ) {
    coupling_factor = kGF2 * P_lep / (32 * kPi * kPi * p_probe);
    fFormFactors.SetModel( fNCFormFactorsModel );
  }
  else if ( proc_info.IsEM() ) {
    coupling_factor = kAem2 * E_lep/ (probeP4.E() * Q2*Q2) ;
    fFormFactors.SetModel( fEMFormFactorsModel );
  }
  else {
    LOG("UnifiedQELPXSec", pERROR) << "Unrecognized process type encountered"
      << " in genie::UnifiedQELPXSec::XSec()";
    return 0.;
  }
 
  // If just testing phase space, set coupling factor to 1
  if(testing_phase_space) coupling_factor = 1.;

  // Apply the coupling factor to the differential cross section
  xsec *= coupling_factor;

  // For a bound hit nucleon, the corrected energy transfer qTilde0 needs to be
  // positive to enforce energy conservation
  if ( qTildeP4.E() <= 0. && interaction->InitState().Tgt().IsNucleus()
    && !interaction->TestBit(kIAssumeFreeNucleon) ) return 0.;

  // Set Q2 to Q2tilde while computing form factors
  interaction->KinePtr()->SetQ2( Q2tilde );

  // Evaluate the form factors
  fFormFactors.Calculate( interaction );

  double f1v = fFormFactors.F1V();
  double f2v = fFormFactors.xiF2V();
  double ffa = fFormFactors.FA();
  double ffp = fFormFactors.Fp();
  
  // Now that we've calculated them, store the true Q2 value
  interaction->KinePtr()->SetQ2( Q2 );

  if(!testing_phase_space) {
    // Compute the tensor contraction
    // Compute the leptonic tensor
    // Note we have to pass a bool to LeptonTensor to use the SF conventions
    InteractionType_t type = interaction->ProcInfo().InteractionTypeId();
    LeptonTensor L_munu( probeP4, lepP4, init_state.ProbePdg(), type, true);

    // Set up the 4x4 hadronic response tensor
    std::complex<double> HadronTensor[4][4]; 

    // xmn_in (GeV)
    double pm_GeV = genie::PDGLibrary::Instance()->Find(2212)->Mass();
    double nm_GeV = genie::PDGLibrary::Instance()->Find(2112)->Mass();
    double xmn_in = 0.5 * (pm_GeV + nm_GeV);

    // Set up dirac matrices routine 
    diracmatrices_(&xmn_in);

    // Get energy and momentum transfer values
    double q = qP4.P();
    double w = qP4.E();
    double wt = qTildeP4.E(); 
    double xk_x, xk_y, xk_z, xp_x, xp_y, xp_z;
    //fill_kin_vars(xk_x, xk_y, xk_z, xp_x, xp_y, xp_z, p4Ni.Vect(), p4Nf.Vect());
    xk_x = p4NiOnShell.X();
    xk_y = p4NiOnShell.Y();
    xk_z = p4NiOnShell.Z();
    xp_x = p4Nf.X();
    xp_y = p4Nf.Y();
    xp_z = p4Nf.Z();  

    // Compute hadron tensor 
    compute_hadron_tensor_(&wt, &xk_x, &xk_y, &xk_z, &xp_x, &xp_y, &xp_z, &f1v, &f2v, &ffa, &ffp, HadronTensor);
  
    // Convert to manual response tensor
    ManualResponseTensor ATilde_munu(HadronTensor);
  /*  
    LOG("UnifiedQELPXsec", pWARN) << "Pi = (" << p4NiOnShell.E() << "," << p4NiOnShell.X() << "," <<p4NiOnShell.Y() << "," <<p4NiOnShell.Z() << ")";
    LOG("UnifiedQELPXsec", pWARN) << "Pf = (" << p4Nf.E() << "," << p4Nf.X() << "," << p4Nf.Y() << "," << p4Nf.Z() << ")";
    LOG("UnifiedQELPXsec", pWARN) << "q = (" << qTildeP4.E() << "," << qTildeP4.X() << "," << qTildeP4.Y() << "," << qTildeP4.Z() << ")";
    LOG("UnifiedQELPXsec", pWARN) << "k = (" << probeP4.E() << "," << probeP4.X() << "," << probeP4.Y() << "," << probeP4.Z() << ")";
    LOG("UnifiedQELPXsec", pWARN) << "kf = (" << lepP4.E() << "," << lepP4.X() << "," << lepP4.Y() << "," << lepP4.Z() << ")";
    LOG("UnifiedQELPXsec", pWARN) << "Q2tilde = " << Q2tilde << ", f1v = " << f1v << ", f2v = " << f2v << ", fa = " << ffa;
    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        LOG("UnifiedQELPXsec", pWARN) << "R(" << i << "," << j << ") = " << ATilde_munu(static_cast<genie::TensorIndex_t>(i),static_cast<genie::TensorIndex_t>(j)) << " ";
      }
    }

    for(int i = 0; i < 4; i++) {
      for(int j = 0; j < 4; j++) {
        LOG("UnifiedQELPXsec", pWARN) << "L(" << i << "," << j << ") = " << L_munu(static_cast<genie::TensorIndex_t>(i),static_cast<genie::TensorIndex_t>(j)) << " ";
      }
    }
  */
    std::complex<double> contraction = L_munu * ATilde_munu;

    //LOG("UnifiedQELPXsec", pWARN) << "Contraction = " << contraction << "\n";

    if ( std::abs(contraction.imag()) > kASmallNum ) {
      LOG("UnifiedQELPXSec", pWARN) << "Tensor contraction has nonvanishing imaginary part!";
    }
  
    // Apply the tensor contraction to the cross section
    xsec *= contraction.real();
  
  }
  // Multiply by the analytic solution of the energy-conserving delta function
  // used by the kPSQELEvGen phase space. 
  xsec *= genie::utils::EnergyDeltaFunctionSolutionQEL( *interaction );
 

  // Check whether variable tranformation is needed
  if ( kps != kPSQELEvGen ) {
    // Compute the appropriate Jacobian for transformation to the requested
    // phase space
    double J = utils::kinematics::Jacobian(interaction, kPSQELEvGen, kps);
    xsec *= J;
  }

  return xsec;
}
//____________________________________________________________________________
double UnifiedQELPXSec::Integral(const Interaction* in) const
{
  // Intended for use with genie::NewQELXSec, which is smart
  // enough to handle free nucleon vs. nuclear targets, different
  // nuclear models, etc.
  return fXSecIntegrator->Integrate(this, in);
}
//____________________________________________________________________________
bool UnifiedQELPXSec::ValidProcess(const Interaction* interaction) const
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

  return false;
}
//____________________________________________________________________________
void UnifiedQELPXSec::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void UnifiedQELPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void UnifiedQELPXSec::LoadConfig(void)
{
  double thc;
  GetParam( "CabibboAngle", thc ) ;
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

  // cross section scaling factor
  GetParam( "QEL-CC-XSecScale", fXSecScale ) ;

  // load QEL form factors models
  fCCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("CCFormFactorsAlg") );
  assert( fCCFormFactorsModel );

  fNCFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("NCFormFactorsAlg") );
  assert( fNCFormFactorsModel );

  fEMFormFactorsModel = dynamic_cast<const QELFormFactorsModelI*>(
    this->SubAlg("EMFormFactorsAlg") );
  assert( fEMFormFactorsModel );

  // Attach CC model for now. This will be updated later.
  fFormFactors.SetModel( fCCFormFactorsModel );

  // load xsec integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI*>(
    this->SubAlg("XSec-Integrator") );
  assert(fXSecIntegrator);

  // get nuclear model
  fNuclModel = dynamic_cast<const NuclearModelI*>(
    this->SubAlg("IntegralNuclearModel") );
  assert(fNuclModel);

  // get the algorithm ID for the PauliBlocker
  RgAlg pauliBlockID;
  GetParamDef( "PauliBlockerAlg", pauliBlockID,
    RgAlg("genie::PauliBlocker", "Default") );
  AlgId pbID = AlgId( pauliBlockID );

  AlgFactory* algf = AlgFactory::Instance();
  fPauliBlocker = dynamic_cast<const PauliBlocker*>(
    algf->GetAlgorithm(pbID) );
  assert( fPauliBlocker );
}
//____________________________________________________________________________
double UnifiedQELPXSec::FreeNucXSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  // This cross section calculation works in the kPSQ2fE phase space
  // and is intended for direct tests against other QEL cross section
  // models using free nucleon targets.

  // Get kinematics and init-state parameters
  const Kinematics&   kinematics = interaction->Kine();
  const InitialState& init_state = interaction->InitState();
  const Target& target = init_state.Tgt();

  // Work in the rest frame of the initial struck nucleon
  TLorentzVector* temp_probeP4 = init_state.GetProbeP4( kRfHitNucRest );
  TLorentzVector probeP4 = *temp_probeP4;
  delete temp_probeP4;
  double E_probe = probeP4.E();
  double p_probe = probeP4.P();
  double m_probe2 = probeP4.M2();

  // Check that we are above threshold for the reaction using Mandelstam s
  double sqrt_s = init_state.CMEnergy();
  double ml = interaction->FSPrimLepton()->Mass();
  double W = interaction->RecoilNucleon()->Mass();
  if ( sqrt_s < W + ml ) return 0.;

  // Hit nucleon mass (may be off-shell)
  double M = target.HitNucP4().M();
  // 4-momentum of initial hit nucleon in its rest frame
  // TODO: think more about bound nucleons here
  TLorentzVector p4Ni(0., 0., 0., M);

  // Use Q^2 to compute the full 4-momenta needed for the tensor contraction
  double Q2 = kinematics.Q2();

  // In the hit nucleon rest frame, we can compute the energy transfer easily
  // from Q^2
  double q0 = ( Q2 + W*W - M*M ) / ( 2. * M );

  // Final lepton total energy and momentum
  double El = E_probe - q0;
  double pl = std::sqrt( std::max(0., El*El - ml*ml) );

  // Also compute the scattering cosine in the hit nucleon rest frame
  // using Q^2
  double cos_theta_l = ( 2.*E_probe*El - Q2 - m_probe2 - ml*ml )
    / ( 2. * p_probe * pl );
  double theta_l = std::acos( cos_theta_l );

  // The cross section is invariant with respect to the azimuthal angle
  // of the outgoing lepton, so just pick zero for simplicity
  double phi_l = 0.;

  // Build the final lepton 4-momentum and 4-momentum transfer
  TLorentzVector p4l(0., 0., pl, El);
  p4l.SetTheta( theta_l );
  p4l.SetPhi( phi_l );

  TLorentzVector qP4 = probeP4 - p4l;

  // Start xsec calculation with overall phase space factor and
  // scaling factor from XML
  double xsec = fXSecScale / ( 64. * kPi * std::pow(M * p_probe, 2) );

  // Get the correct couplings and form factors model for the current
  // interaction
  // TODO: reduce code duplication here
  double coupling_factor = 1.;
  const genie::ProcessInfo& proc_info = interaction->ProcInfo();
  if ( proc_info.IsWeakCC() ) {
    coupling_factor = kGF2 * fCos8c2 / 2.;
    fFormFactors.SetModel( fCCFormFactorsModel );
  }
  else if ( proc_info.IsWeakNC() ) {
    coupling_factor = kGF2 / 2.;
    fFormFactors.SetModel( fNCFormFactorsModel );
  }
  else if ( proc_info.IsEM() ) {
    // Lorentz-Heaviside natural units: e^4 = 16 * kPi^2 * kAem2
    coupling_factor = 16. * kPi * kPi * kAem2 / ( Q2*Q2 );
    fFormFactors.SetModel( fEMFormFactorsModel );
  }
  else {
    LOG("UnifiedQELPXSec", pERROR) << "Unrecognized process type encountered"
      << " in genie::CBFSpecFuncQELPXSec::XSec()";
    return 0.;
  }

  // Apply the coupling factor to the differential cross section
  xsec *= coupling_factor;

  // Evaluate the form factors
  fFormFactors.Calculate( interaction );

  // Compute the tensor contraction
  InteractionType_t type = interaction->ProcInfo().InteractionTypeId();
  int hit_nuc_pdg = target.HitNucPdg();

  LeptonTensor L_munu( probeP4, p4l, init_state.ProbePdg(), type );

  // TO DO:: Implement free nucleon tensor with fotran wrapper
  std::complex<double> HadronTensor[4][4];
  // Convert to manual response tensor
  ManualResponseTensor ATilde_munu(HadronTensor);

  std::complex<double> contraction = L_munu * ATilde_munu;

  if ( std::abs(contraction.imag()) > kASmallNum ) {
    LOG("UnifiedQELPXSec", pWARN) << "Tensor contraction has nonvanishing"
      << " imaginary part!";
  }

  // Apply the tensor contraction to the cross section
  xsec *= contraction.real();

  // Check whether variable transformation is needed
  if ( kps != kPSQ2fE ) {
    // Compute the appropriate Jacobian for transformation to the requested
    // phase space
    double J = utils::kinematics::Jacobian(interaction, kPSQ2fE, kps);
    xsec *= J;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt
  if ( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear suppression factor
  // (R(Q2) is adapted from NeuGEN - see comments therein)
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);

  // Number of scattering centers in the target
  int NNucl = pdg::IsProton(hit_nuc_pdg) ? target.Z() : target.N();

  LOG("UnifiedQELPXSec", pDEBUG) << "Nuclear suppression factor R(Q2) = "
    << R << ", NNucl = " << NNucl;

  // Compute nuclear cross section
  xsec *= ( R * NNucl );

  return xsec;
}

