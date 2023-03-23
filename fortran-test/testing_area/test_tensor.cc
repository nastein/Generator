// Run with genie -l to load the needed class dictionaries into ROOT
#include <fstream>
#include <iostream>
//#include <iomanip>

#include "Physics/QuasiElastic/XSection/ManualResponseTensor.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/QuasiElastic/XSection/FreeNucleonTensor.h"
#include "Physics/QuasiElastic/XSection/LeptonTensor.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/SpectralFunc.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/RunOpt.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Conventions/RefFrame.h"


//using namespace std;

//extern"C" {
//void fortfunc_(std::complex<double> array[4][4]);
//}

//int main(int argc, char *argv[]) {
int test_tensor() {
/*
  std::complex<double> tensor[4][4];
  cout << "Initial array\n";
  for(int i = 0; i < 4; i++) {
      cout << array[i][0] << ", " << array[i][1] << ", " << array[i][2] << array[i][3] <<"\n";
   }
  
  fortfunc_(tensor);
  cout << "Final array\n";
  for(int i = 0; i < 4; i++) {
      cout << array[i][0] << ", " << array[i][1] << ", " << array[i][2] << array[i][3] <<"\n";
   }
*/

  std::complex<double> tensor[4][4] = {{2,0,0,0},{0,2,0,0},{0,0,2,0},{0,0,0,2}};

  genie::ManualResponseTensor W_munu(tensor);
  std::cout << "Manual response tensor array\n" << std::endl;
  for(int i = 0; i < 4; i++) {
      std::cout <<  W_munu(static_cast<genie::TensorIndex_t>( i ), static_cast<genie::TensorIndex_t>( 0 )) << ", " <<  W_munu(static_cast<genie::TensorIndex_t>( i ), static_cast<genie::TensorIndex_t>( 1 )) << ", " <<  W_munu(static_cast<genie::TensorIndex_t>( i ), static_cast<genie::TensorIndex_t>( 2 )) << ", " <<  W_munu(static_cast<genie::TensorIndex_t>( i ), static_cast<genie::TensorIndex_t>( 3 )) <<"\n" << std::endl;
   };

  std::cout << "Contracted tensor" << std::endl;
  std::cout << W_munu * W_munu << std::endl;

  return 0;

}
