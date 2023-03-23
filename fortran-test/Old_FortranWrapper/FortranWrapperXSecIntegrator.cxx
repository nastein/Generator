//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/Common/VertexGenerator.h"
#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelMap.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TGenPhaseSpace.h"
#include "FortranWrapperXSecIntegrator.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::gsl;

//____________________________________________________________________________
FortranWrapperXSecIntegrator::FortranWrapperXSecIntegrator() : XSecIntegratorI("genie::FortranWrapperXSecIntegrator")
{

}
//____________________________________________________________________________
FortranWrapperXSecIntegrator::FortranWrapperXSecIntegrator(std::string config) : XSecIntegratorI("genie::FortranWrapperXSecIntegrator", config)
{

}
//____________________________________________________________________________
double FortranWrapperXSecIntegrator::Integrate(const XSecAlgorithmI* model, const Interaction* in) const
{
  LOG("FortranWrapperXSecIntegrator",pDEBUG) << "Beginning integrate";
  if ( !model->ValidProcess(in) ) return 0.;

  Interaction* interaction = new Interaction( *in );
  interaction->SetBit( kISkipProcessChk );
  interaction->SetBit( kISkipKinematicChk );

  const NuclearModelI* nucl_model = dynamic_cast< const NuclearModelI* >(
    model->SubAlg("NuclearModel") );
  assert( nucl_model );

  AlgFactory* algf = AlgFactory::Instance();
  const VertexGenerator* vtx_gen = dynamic_cast<const VertexGenerator*>(
    algf->GetAlgorithm(fVertexGenID) );
  assert( vtx_gen );

  // Determine the appropriate binding energy mode to use.
  // The default given here is for the case of a free nucleon.
  QELEvGen_BindingMode_t bind_mode = kOnShell;
  Target* tgt = interaction->InitState().TgtPtr();
  if ( tgt->IsNucleus() ) {
    std::string bind_mode_str = model->GetConfig()
      .GetString( "IntegralNucleonBindingMode" );
    bind_mode = genie::utils::StringToQELBindingMode( bind_mode_str );
  }

  // Set the masses of the two-body final state for QE scattering. These are
  // needed for phase-space generation below
  const double ml = interaction->FSPrimLepton()->Mass();
  const double mNf = interaction->RecoilNucleon()->Mass();
  double masses[2] = { ml, mNf };
  //std::cout << "m[0] = " << ml << "\n";
  //std::cout << "m[1] = " << mNf << "\n";
  

  double xsec_max = 0.;
  double xsec_sum_of_squares = 0.;
  double xsec_sum = 0.;
  bool converged = false;
  unsigned int n = 0u;

  fMaxDiffXSec = -1.;
  while ( !converged && n < fMaxThrows ) {

    ++n;

    // For a free nucleon target (hit nucleon is at rest in the lab frame), we
    // don't need to do an MC integration over the initial state variables. In
    // this case, just set up the nucleon at the origin, on-shell, and at rest.
    if ( !tgt->IsNucleus() ) {
      tgt->SetHitNucPosition(0.);
      nucl_model->SetMomentum3( TVector3(0., 0., 0.) );
      nucl_model->SetRemovalEnergy(0.);
    }
    // For a nuclear target, we need to loop over a bunch of nucleons sampled
    // from the nuclear model (with positions sampled from the vertex generator
    // to allow for using the local Fermi gas model).
    else {
      // Select a new position for the initial hit nucleon (needed for the local
      // Fermi gas model, but other than slowing things down a bit, it doesn't
      // hurt to do this for other models)
      TVector3 vertex_pos = vtx_gen->GenerateVertex( interaction, tgt->A() );
      double radius = vertex_pos.Mag();
      tgt->SetHitNucPosition( radius );

      // Sample a new nucleon 3-momentum and removal energy (this will be applied
      // to the nucleon via a call to genie::utils::ComputeFullQELPXSec(), so
      // there's no need to mess with its 4-momentum here)
      nucl_model->GenerateNucleon(*tgt, radius);
    }

    // Set the initial struck nucleon 4-momentum appropriately based on the
    // information from the nuclear model
    double dummy_Ebind;
    genie::utils::BindHitNucleon( *interaction, *nucl_model,
      dummy_Ebind, bind_mode );

    // TODO: reduce code duplication between this function and FortranWrapperEventGenerator
    // Check that we're above threshold for the QE reaction. If we're not, then
    // start over with a new nucleon
    double sqrt_s = interaction->InitState().CMEnergy();
    double Ediff = sqrt_s - ml - mNf;
    if ( Ediff < 0. || sqrt_s < 0. || isnan(sqrt_s) || isnan(Ediff) ) continue;

    // Determine the total 4-momentum of the two-body system (probe, bound
    // nucleon) in the initial state
    TLorentzVector p4Ni = interaction->InitState().Tgt().HitNucP4();
    TLorentzVector* p4ProbeTemp = interaction->InitState().GetProbeP4( kRfLab );
    TLorentzVector p4Probe = *p4ProbeTemp;
    delete p4ProbeTemp;
    TLorentzVector p4tot = p4Probe + p4Ni;
  
    // Set up the phase-space generator
    RandomGen * rnd = RandomGen::Instance();
    TGenPhaseSpace ps_gen;
    bool permitted = ps_gen.SetDecay( p4tot, 2, masses );
    if(!permitted) {
	std::cout << "Phase space decay not permitted!" << "\n";
	double sqrt_s = interaction->InitState().CMEnergy();
	std::cout << "sqrt_s = " << sqrt_s << "\n";   
 	double Ediff = sqrt_s - ml - mNf;
    	std::cout << "Ediff = " << Ediff << "\n";
	std::cout << "CM energy = " << p4tot.Mag() << "\n";
	std::cout << "P4Probe = " << "\n";
	p4Probe.Print();
        std::cout << "P4Ni = " <<  "\n";
	p4Ni.Print();
        std::cout << "P4tot = " << "\n";
	p4tot.Print();
	std::cout << "ml = " << ml << "\n";
	std::cout << "mNf = " << mNf << "\n";
        std::cout << "CM energy - sum of masses = " << p4tot.Mag() - ml - mNf << "\n";
	throw; 
    }

    double weight_max_ps = ps_gen.GetWtMax();
    //std::cout << "Max ps weight = " << weight_max_ps << "\n";
    bool accept_decay=false;
    while(!accept_decay) 
    {
    	// Choose a phase space point uniformly using the generator. We don't need
    	// to use accept/reject here because a suitable point is always found for a
    	// two-body system.
    	double weight = ps_gen.Generate();
    	double gw = weight_max_ps * rnd->RndDec().Rndm();
    	accept_decay = (gw<=weight);
     }
         
   


    // Retrieve the final particle 4-momenta from the phase-space generator
    TLorentzVector* p4Lep = ps_gen.GetDecay( 0 );
    TLorentzVector* p4Nf  = ps_gen.GetDecay( 1 );

    // Update the interaction object with the sampled values
    interaction->KinePtr()->SetFSLeptonP4( *p4Lep );
    interaction->KinePtr()->SetHadSystP4( *p4Nf );

    // Compute the differential cross section
    double xsec = model->XSec( interaction, kPSFullNBody);

    if ( std::isnan(xsec) ) continue;

    if ( xsec > fMaxDiffXSec ) fMaxDiffXSec = xsec;

    //xsec *= (kPi/(2.*pow(sqrt_s*2.*kPi,2)))*sqrt(()*());  

    xsec_sum += xsec;

    xsec_sum_of_squares += xsec * xsec;

    if ( xsec_sum > 0. && n > 2000000 ) {
      double xsec_mean = xsec_sum / n;
      double sum_variance = (xsec_sum_of_squares / n) - xsec_mean*xsec_mean;
      double xsec_std_error = std::sqrt( sum_variance / n );
      double rel_error = xsec_std_error / xsec_mean;
      if ( n % 200000 == 0 ) {
        std::cout << "  ITERATION = " << n << '\n';
        std::cout << "  MAX DIFF XSEC = " << fMaxDiffXSec << '\n';
        std::cout << "  INTEGRAL = " << xsec_mean
          << " +/- " << xsec_std_error << '\n';
        std::cout << "  SUM VARIANCE = " << sum_variance << '\n';
        std::cout << "  STD ERROR = " << xsec_std_error << '\n';
        std::cout << "  REL ERROR = " << rel_error << "\n\n\n";
      }
      if ( rel_error < fDesiredRelativePrecision ) converged = true;
    }
  }

  // Delete our temporary clone of the input Interaction
  delete interaction;

  // MC estimator of the total cross section is the mean of the xsec values
  double xsec_mean = xsec_sum / n;

  std::cout << "Integral = " << xsec_mean << ", iterations = " << n << '\n';

  return xsec_mean;
}
//____________________________________________________________________________
void FortranWrapperXSecIntegrator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperXSecIntegrator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void FortranWrapperXSecIntegrator::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  int max;
  GetParamDef( "MaxThrows", max, 5000000 ) ;
  fMaxThrows = static_cast<unsigned int>( max );

  GetParamDef( "DesiredRelativePrecision", fDesiredRelativePrecision, 0.001 );

  RgAlg vertexGenID;
  GetParamDef( "VertexGenAlg", vertexGenID, RgAlg("genie::VertexGenerator", "Default") );
  fVertexGenID = AlgId( vertexGenID );

  fMaxDiffXSec = 0.;
}
