// Run with genie -l to load the needed class dictionaries into ROOT
#include <fstream>
#include <iostream>
//#include <iomanip>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Interaction/InteractionType.h"
#include "Physics/Fortran/XSection/UnifiedQELPXSec.h"

//using namespace genie;

constexpr int PROTON = 2212; // proton
constexpr int NEUTRON = 2112; // neutron
constexpr long int rng_seed = 41285;

double testfort(double Ev, long int seed) {

  // Initialize GENIE with the appropriate tune
  genie::RunOpt* run_opt = genie::RunOpt::Instance();

  run_opt->SetTuneName( "N19_00a_00_000" );
  run_opt->BuildTune();
  std::cout << "Random number = " << seed << "\n";
  genie::utils::app_init::RandGen( seed );

  int target_pdg = 1000060120; // 40Ar
  int probe_pdg = 14; // nu_mu

  double proton_xsec = 0.;
  double neutron_xsec = 0.;

  genie::AlgFactory* factory = genie::AlgFactory::Instance();

  const auto* xsec_model = dynamic_cast<const genie::XSecAlgorithmI*>(
    factory->GetAlgorithm("genie::UnifiedQELPXSec", "Dipole") );
  
  //genie::Interaction* interaction_p = genie::Interaction::QELCC( target_pdg, PROTON,
   // probe_pdg, Ev );
  genie::Interaction* interaction_n = genie::Interaction::QELCC( target_pdg, NEUTRON,
    probe_pdg, Ev );

  //proton_xsec = xsec_model->Integral( interaction_p );
  neutron_xsec = xsec_model->Integral( interaction_n );

  //std::cout << "Proton Final result: E = "<< Ev  << ", xsec = " << proton_xsec << " GeV^(-2)\n";
  std::cout << "Neutron Final result: E = "<< Ev  << ", xsec = " << neutron_xsec << " GeV^(-2)\n";
  // Convert from GeV^(-2) to nb
  neutron_xsec /= genie::units::nanobarn;

  std::cout << neutron_xsec << " nb = " << "\n";

  // Convert from nb to 10^(-38) cm^2
  neutron_xsec *=1e5;

  std::cout << neutron_xsec << " 1e-38 cm^(2)\n";

  return 0;
}

