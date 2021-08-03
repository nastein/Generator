//____________________________________________________________________________
/*!

\class    genie::MECXSec

\brief    A numerical cross-section integrator (GENIE/GSL interface) for the
          J. Nieves, I. Ruiz Simo, M.J. Vicente Vacas and Martini MEC models.
          Is a concrete implementation of the XSecIntegratorI interface.

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  March 22, 2016

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _MEC_XSEC_H_
#define _MEC_XSEC_H_

#include "Physics/XSectionIntegration/XSecIntegratorI.h"

#include <Math/Integrator.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>

namespace genie {

class XSecAlgorithmI;
class Interaction;

class MECXSec : public XSecIntegratorI {
public:
  MECXSec();
  MECXSec(string config);
  virtual ~MECXSec();

  // XSecIntegratorI interface implementation
  double Integrate(const XSecAlgorithmI * model, const Interaction * i) const;

  // Overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

protected:
  bool fSplitIntegral;

private:
  void LoadConfig (void);
  double fQ3Max;
};

} // genie namespace

#endif  // _MEC_XSEC_H_
