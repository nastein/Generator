//____________________________________________________________________________
/*!

\class    genie::FortranWrapperEventGenerator

\brief    Generates values for the kinematic variables describing QEL neutrino
          interaction events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Steven Gardiner
          Fermi National Accelerator Laboratory <gardiner \at fnal.gov>

\created  July 26, 2021

\cpright  Copyright (c) 2003-2021, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _FORTRAN_WRAPPER_GENERATOR_H_
#define _FORTRAN_WRAPPER_GENERATOR_H_

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/Common/KineGeneratorWithCache.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Conventions/Controls.h"

namespace genie {

class FortranWrapperEventGenerator: public KineGeneratorWithCache {

public :
  FortranWrapperEventGenerator();
  FortranWrapperEventGenerator(std::string config);
 ~FortranWrapperEventGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord* event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry& config);
  void Configure(std::string config);

private:

  void   LoadConfig     (void);
  double ComputeMaxXSec(const Interaction* in) const;

  /// Add the recoiling spectator nucleus to the event record
  void AddTargetNucleusRemnant(GHepRecord* evrec) const;

  /// Nuclear model used to simulate the initial struck nucleon
  const NuclearModelI* fNuclModel;

  /// Enum that indicates which approach should be used to handle the binding
  /// energy of the struck nucleon
  QELEvGen_BindingMode_t fHitNucleonBindingMode;

}; // class definition

} // genie namespace

#endif // _FORTRAN_WRAPPER_GENERATOR_H_
