//____________________________________________________________________________
/*!

\class    genie::FortranWrapperXSecIntegrator

\brief    Computes the Quasi Elastic (QEL) total cross section. \n
          Is a concrete implementation of the XSecIntegratorI interface. \n

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  July 26, 2021

\cpright  Copyright (c) 2003-2021, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _FORTRAN_WRAPPER_XSEC_INTEGRATOR_H_
#define _FORTRAN_WRAPPER_XSEC_INTEGRATOR_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELUtils.h"

#include "TMath.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

namespace genie {

class NuclearModelI;
class VertexGenerator;

class FortranWrapperXSecIntegrator : public XSecIntegratorI {

public:

  FortranWrapperXSecIntegrator(void);
  FortranWrapperXSecIntegrator(std::string config);

  /// XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI* model, const Interaction* i) const;

  /// Overload the Algorithm::Configure() methods to load private data
  /// members from configuration options
  void Configure(const Registry& config);
  void Configure(std::string config);

  inline double GetMaxDiffXSec() const { return fMaxDiffXSec; }

private:

  void LoadConfig (void);

  // XML configuration parameters
  unsigned int fMaxThrows;
  double fDesiredRelativePrecision;
  AlgId fVertexGenID;
  mutable double fMaxDiffXSec;
};


} // genie namespace

#endif  // _FORTRAN_WRAPPER_XSEC_INTEGRATOR_H_
