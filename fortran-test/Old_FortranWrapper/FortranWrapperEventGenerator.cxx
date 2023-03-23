//____________________________________________________________________________
/*
 Copyright (c) 2003-2021, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/Common/PrimaryLeptonUtils.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/PrintUtils.h"

// Extra includes
#include "TGenPhaseSpace.h"
#include "FortranWrapperEventGenerator.h"
#include "FortranWrapperXSecIntegrator.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
FortranWrapperEventGenerator::FortranWrapperEventGenerator() :
  KineGeneratorWithCache("genie::FortranWrapperEventGenerator")
{

}
//___________________________________________________________________________
FortranWrapperEventGenerator::FortranWrapperEventGenerator(string config) :
  KineGeneratorWithCache("genie::FortranWrapperEventGenerator", config)
{

}
//___________________________________________________________________________
FortranWrapperEventGenerator::~FortranWrapperEventGenerator()
{

}
//___________________________________________________________________________
void FortranWrapperEventGenerator::ProcessEventRecord(GHepRecord* evrec) const
{
  LOG("QELEvent", pDEBUG) << "Generating QE event kinematics...";
  // Get access to the random number generator
  RandomGen* rnd = RandomGen::Instance();

  // Access cross section algorithm for running thread
  RunningThreadInfo* rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI* evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  // Get the interaction and check that we are working with a nuclear target
  Interaction* interaction = evrec->Summary();
  // Skip if not a nuclear target
  if ( interaction->InitState().Tgt().IsNucleus() ) {
    // Skip if no hit nucleon is set
    if ( !evrec->HitNucleon() ) {
      LOG("QELEvent", pFATAL) << "No hit nucleon was set";
      gAbortingInErr = true;
      std::exit(1);
    }
  } // is nuclear target

  // set the 'trust' bits
  interaction->SetBit( kISkipProcessChk );
  interaction->SetBit( kISkipKinematicChk );

  // Access the hit nucleon and target nucleus entries at the GHEP record
  GHepParticle* nucleon = evrec->HitNucleon();
  GHepParticle* nucleus = evrec->TargetNucleus();
  bool have_nucleus = ( nucleus != 0 );

  // Store the hit nucleon radius before computing the maximum differential
  // cross section (important when using the local Fermi gas model)
  Target* tgt = interaction->InitState().TgtPtr();
  double hitNucPos = nucleon->X4()->Vect().Mag();
  tgt->SetHitNucPosition( hitNucPos );

  // Get the maximum differential cross section to use in rejection sampling
  // below
  double xsec_max = this->MaxXSec( evrec );

  // For a composite nuclear target, check to make sure that the
  // final nucleus has a recognized PDG code
  if ( have_nucleus ) {
    // compute A,Z for final state nucleus & get its PDG code
    int nucleon_pdgc = nucleon->Pdg();
    bool is_p  = pdg::IsProton( nucleon_pdgc );
    int Zi = nucleus->Z();
    int Z = ( is_p ) ? Zi - 1 : Zi;
    int A = nucleus->A() - 1;
    TParticlePDG* fnucleus = 0;
    int ipdgc = pdg::IonPdgCode( A, Z );
    fnucleus = PDGLibrary::Instance()->Find( ipdgc );
    if ( !fnucleus ) {
      LOG("QELEvent", pFATAL) << "No particle with [A = " << A << ", Z = " << Z
        << ", pdgc = " << ipdgc << "] in PDGLibrary!";
      std::exit(1);
    }
  }

  // Set the masses of the two-body final state for QE scattering. These are
  // needed for phase-space generation below
  const double ml = interaction->FSPrimLepton()->Mass();
  const double mNf = interaction->RecoilNucleon()->Mass();
  double masses[2] = { ml, mNf };

  // In the accept/reject loop, each iteration samples new values of the
  // initial struck nucleon 4-momentum and the 4-momenta of the outgoing lepton
  // and nucleon. All of these are evaluated in the laboratory frame.
  unsigned int iter = 0;
  bool accept = false;
  while ( true ) {

    ++iter;
    LOG("QELEvent", pINFO) << "Attempt #: " << iter;
    if ( iter > kRjMaxIterations ) {
      LOG("QELEvent", pWARN)
        << "Couldn't select valid event kinematics after "
        << iter << " iterations";
      evrec->EventFlags()->SetBitNumber( kKineGenErr, true );
      genie::exceptions::EVGThreadException exception;
      exception.SetReason( "Couldn't select kinematics" );
      exception.SwitchOnFastForward();
      throw exception;
    }

    // If the target is a composite nucleus, then sample an initial nucleon
    // 3-momentum and removal energy from the nuclear model.
    
    if ( tgt->IsNucleus() ) {
      fNuclModel->GenerateNucleon( *tgt, hitNucPos );
    }
    
    //Commented out to do freenucleon cross sections 
    else {
      // Otherwise, just set the nucleon to be at rest in the lab frame and
      // unbound. Use the nuclear model to make these assignments. The call
      // to BindHitNucleon() will apply them correctly below.
      fNuclModel->SetMomentum3( TVector3(0., 0., 0.) );
      fNuclModel->SetRemovalEnergy( 0. );
    }

    // Set the initial struck nucleon 4-momentum appropriately based on the
    // information from the nuclear model
    double fEb;
    genie::utils::BindHitNucleon( *interaction, *fNuclModel,
      fEb, fHitNucleonBindingMode);

    // Check that we're above threshold for the QE reaction. If we're not, then
    // start over with a new nucleon
    double sqrt_s = interaction->InitState().CMEnergy();
    double Ediff = sqrt_s - ml - mNf;
    if ( Ediff < 0. || isnan(Ediff) || isnan(sqrt_s)) continue;

    // Determine the total 4-momentum of the two-body system (probe, bound
    // nucleon) in the initial state
    TLorentzVector p4Ni = interaction->InitState().Tgt().HitNucP4();
    TLorentzVector* p4ProbeTemp = interaction->InitState().GetProbeP4( kRfLab );
    TLorentzVector p4Probe = *p4ProbeTemp;
    delete p4ProbeTemp;
    TLorentzVector p4tot = p4Probe + p4Ni;

    // Set up the phase-space generator
    TGenPhaseSpace ps_gen;
    ps_gen.SetDecay( p4tot, 2, masses );

    // Maximum weight for sampling over the allowed phase space
    double max_weight_ps = ps_gen.GetWtMax();

    // Choose a phase space point uniformly using the generator. We don't need
    // to use accept/reject here because a suitable point is always found for a
    // two-body system.
    double weight = ps_gen.Generate();
    //std::cout << "weight = " << weight << "\n";
    // Retrieve the final particle 4-momenta from the phase-space generator
    TLorentzVector* p4Lep = ps_gen.GetDecay( 0 );
    TLorentzVector* p4Nf  = ps_gen.GetDecay( 1 );

    // Update the interaction object with the sampled values
    interaction->KinePtr()->SetFSLeptonP4( *p4Lep );
    interaction->KinePtr()->SetHadSystP4( *p4Nf );

    // Compute the differential cross section
    double xsec = fXSecModel->XSec( interaction, kPSFullNBody );

   // xsec *= 1. / ( 4. * p4Lep->E() * p4Nf->E() );

    // Accept/reject the current event
    this->AssertXSecLimits( interaction, xsec, xsec_max );

    double t = xsec_max * rnd->RndKine().Rndm();

    accept = ( t < xsec );

    // If the generated kinematics are accepted, finish-up module's job
    if ( accept ) {
      double gQ2 = interaction->KinePtr()->Q2( false );
      LOG("QELEvent", pINFO) << "*Selected* Q^2 = " << gQ2 << " GeV^2";

      // reset bits
      interaction->ResetBit( kISkipProcessChk );
      interaction->ResetBit( kISkipKinematicChk );

      // get neutrino energy at struck nucleon rest frame and the
      // struck nucleon mass (can be off the mass shell)
      const InitialState & init_state = interaction->InitState();
      double E = init_state.ProbeE( kRfHitNucRest );
      double M = init_state.Tgt().HitNucP4().M();
      LOG("QELKinematics", pNOTICE) << "E = " << E << ", M = "<< M;

      // The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
      int rpdgc = interaction->RecoilNucleonPdg();
      double gW = PDGLibrary::Instance()->Find( rpdgc )->Mass();
      LOG("QELEvent", pNOTICE) << "Selected: W = " << gW;

      // (W,Q2) -> (x,y)
      double gx = 0;
      double gy = 0;
      kinematics::WQ2toXY( E, M, gW, gQ2, gx, gy );

      // Lock the selected kinematics & clear running values
      interaction->KinePtr()->SetQ2( gQ2, true );
      interaction->KinePtr()->SetW ( gW,  true );
      interaction->KinePtr()->Setx ( gx,  true );
      interaction->KinePtr()->Sety ( gy,  true );
      interaction->KinePtr()->ClearRunningValues();

      // Set the cross section for the selected kinematics
      evrec->SetDiffXSec( xsec, kPSFullNBody );

      TLorentzVector lepton( interaction->KinePtr()->FSLeptonP4() );
      TLorentzVector outNucleon( interaction->KinePtr()->HadSystP4() );
      TLorentzVector x4l( *(evrec->Probe())->X4() );

      // Add the final-state lepton to the event record
      evrec->AddParticle( interaction->FSPrimLeptonPdg(), kIStStableFinalState,
        evrec->ProbePosition(), -1, -1, -1,
        interaction->KinePtr()->FSLeptonP4(), x4l );

      // Set its polarization
      utils::SetPrimaryLeptonPolarization( evrec );

      // Add the final-state nucleon to the event record
      GHepStatus_t ist = ( tgt->IsNucleus() )
        ? kIStHadronInTheNucleus : kIStStableFinalState;
      evrec->AddParticle( interaction->RecoilNucleonPdg(), ist,
        evrec->HitNucleonPosition(), -1, -1, -1,
        interaction->KinePtr()->HadSystP4(), x4l );

      // Store struck nucleon momentum and binding energy
      TLorentzVector p4ptr = interaction->InitStatePtr()->TgtPtr()->HitNucP4();
      nucleon->SetMomentum( p4ptr );
      nucleon->SetRemovalEnergy( fEb );

      // Add the recoiling remnant nucleus
      this->AddTargetNucleusRemnant( evrec );

      // We've accepted and processed the event, so exit the accept/reject loop
      break;
    }
    else {
      LOG("QELEvent", pDEBUG) << "Reject current throw...";
    }

  } // Accept/reject iterations
  LOG("QELEvent", pINFO) << "Done generating QE event kinematics!";
}
//___________________________________________________________________________
void FortranWrapperEventGenerator::AddTargetNucleusRemnant( GHepRecord* evrec ) const
{
  // add the remnant nuclear target at the GHEP record

  LOG("QELEvent", pINFO) << "Adding final state nucleus";

  double Px = 0;
  double Py = 0;
  double Pz = 0;
  double E  = 0;

  GHepParticle * nucleus = evrec->TargetNucleus();
  bool have_nucleus = nucleus != 0;
  if (!have_nucleus) return;

  int A = nucleus->A();
  int Z = nucleus->Z();

  int fd = nucleus->FirstDaughter();
  int ld = nucleus->LastDaughter();

  for(int id = fd; id <= ld; id++) {

    // compute A,Z for final state nucleus & get its PDG code and its mass
    GHepParticle * particle = evrec->Particle(id);
    assert(particle);
    int  pdgc = particle->Pdg();
    bool is_p  = pdg::IsProton (pdgc);
    bool is_n  = pdg::IsNeutron(pdgc);

    if (is_p) Z--;
    if (is_p || is_n) A--;

    Px += particle->Px();
    Py += particle->Py();
    Pz += particle->Pz();
    E  += particle->E();

  }//daughters

  TParticlePDG * remn = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  remn = PDGLibrary::Instance()->Find(ipdgc);
  if(!remn) {
    LOG("HadronicVtx", pFATAL)
      << "No particle with [A = " << A << ", Z = " << Z
      << ", pdgc = " << ipdgc << "] in PDGLibrary!";
    assert(remn);
  }

  double Mi = nucleus->Mass();
  Px *= -1;
  Py *= -1;
  Pz *= -1;
  E = Mi-E;

  // Add the nucleus to the event record
  LOG("QELEvent", pINFO)
    << "Adding nucleus [A = " << A << ", Z = " << Z
    << ", pdgc = " << ipdgc << "]";

  int imom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
      ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);

  LOG("QELEvent", pINFO) << "Done";
  LOG("QELEvent", pINFO) << *evrec;
}
//___________________________________________________________________________
void FortranWrapperEventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperEventGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperEventGenerator::LoadConfig(void)
{
  // Load sub-algorithms and config data in advance to reduce the number of
  // registry lookups
  fNuclModel = 0;

  RgKey nuclkey = "NuclearModel";

  fNuclModel = dynamic_cast< const NuclearModelI* >(
    this->SubAlg("NuclearModel") );
  assert( fNuclModel );

  // Safety factor for the maximum differential cross section
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.6 );

  // Minimum energy for which max xsec would be cached, forcing explicit
  // calculation for lower eneries
  GetParamDef( "Cache-MinEnergy", fEMin, 1.00 ) ;

  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. );
  assert( fMaxXSecDiffTolerance >= 0. );

  // Decide how to handle the binding energy of the initial state struck
  // nucleon
  std::string binding_mode;
  GetParamDef( "HitNucleonBindingMode", binding_mode, std::string("UseNuclearModel") );

  fHitNucleonBindingMode = genie::utils::StringToQELBindingMode( binding_mode );
}
//____________________________________________________________________________
double FortranWrapperEventGenerator::ComputeMaxXSec(const Interaction* in) const
{
  // Computes the maximum differential cross section in the requested phase
  // space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
  // method and the value is cached at a circular cache branch for retrieval
  // during subsequent event generation.
  // The computed max differential cross section does not need to be the exact
  // maximum. The number used in the rejection method will be scaled up by a
  // safety factor. But it needs to be fast.
  LOG("QELEvent", pINFO) << "Computing maximum cross section to throw against";

  // Use a brute-force scan by calling the integrator and then retrieving the
  // largest differential cross section value seen during MC integration throws
  const FortranWrapperXSecIntegrator* integrator
    = dynamic_cast< const FortranWrapperXSecIntegrator* >(
    fXSecModel->SubAlg( "XSec-Integrator" )
  );
  assert( integrator );

  integrator->Integrate( fXSecModel, in );
  double max = integrator->GetMaxDiffXSec();
  max *= fSafetyFactor;
  return max;
}
