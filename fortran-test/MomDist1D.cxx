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
#include "UnifiedQELPXSec.h"
//using namespace genie;

constexpr int PROTON = 2212; // proton
constexpr int NEUTRON = 2112; // neutron
constexpr long int rng_seed = 41285;
int main(int argc, char *argv[]) {

  // Initialize GENIE with the appropriate tune
  genie::RunOpt* run_opt = genie::RunOpt::Instance();

  TH2D* sf_hist = nullptr;

  run_opt->SetTuneName( "N19_00a_00_000" );
  run_opt->BuildTune();
  std::cout << "Random number = " << strtol(argv[2], NULL, 0) << "\n";
  genie::utils::app_init::RandGen( strtol(argv[2], NULL, 0) );
  
  double Ev = atof(argv[1]); // GeV
  int target_pdg = 1000060120; // 12C
  int probe_pdg = 11; // electron

  double proton_xsec = 0.;
  double neutron_xsec = 0.;

  genie::AlgFactory* factory = genie::AlgFactory::Instance();

  const auto* nuclear_model = dynamic_cast<const genie::SpectralFunc*>(
    factory->GetAlgorithm( "genie::SpectralFunc", "Default" ) );

  genie::Target tmp_tgt( target_pdg );
  sf_hist = nuclear_model->SelectSpectralFunction( tmp_tgt );
  
  TCanvas *c1 = new TCanvas("","",1200,800);
  sf_hist->Draw();
  c1->SaveAs("SF_2D.png");
  
  return 0;
}

