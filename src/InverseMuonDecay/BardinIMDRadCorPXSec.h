//____________________________________________________________________________
/*!

\class    genie::BardinIMDRadCorPXSec

\brief    Computes the Inverse Muon Decay (IMD) diff. cross section using the 
          Bardin-Dokuchaeva including all 1-loop radiative corrections. \n

          This is a 'trully' inclusive IMD cross section, i.e. the brem. cross
          section (dxsec_brem/dy)|w>w0 [see Bardin paper, cited below] is not
          subtracted from the IMD cross section and therefore it is not suitable
          for experimental situations where a photon energy trigger threshold
          is applied.

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      D.Yu.Bardin and V.A.Dokuchaeva, Nucl.Phys.B287:839 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 14, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _BARDIN_IMD_RADIATIVE_CORRECTIONS_PARTIAL_XSEC_H_
#define _BARDIN_IMD_RADIATIVE_CORRECTIONS_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Numerical/GSFunc.h"

namespace genie {

class IntegratorI;

class BardinIMDRadCorPXSec : public XSecAlgorithmI {

public:
  BardinIMDRadCorPXSec();
  BardinIMDRadCorPXSec(string config);
  virtual ~BardinIMDRadCorPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  bool   ValidProcess    (const Interaction * i) const;
  bool   ValidKinematics (const Interaction * i) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  //-- load configuration when Algorithm::Configure() is called
  void LoadConfig(void);

  //-- Private functions
  //   (symbols follow the notation in Bardin-Dokuchaeva paper)
  double Li2 (double z)                      const;
  double Fa  (double re, double r, double y) const;
  double P   (int    i,  double r, double y) const;
  double C   (int    i,  int k,    double r) const;

  //-- data members
  const IntegratorI * fIntegrator;
};

} // genie namespace

//____________________________________________________________________________
/*!
\class    genie::BardinIMDRadCorIntegrand

\brief    Auxiliary scalar function for the internal integration in Bardin's
          IMD d2xsec/dxdy cross section algorithm

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 20, 2006
*/
//____________________________________________________________________________

namespace genie {

class BardinIMDRadCorIntegrand : public GSFunc
{
public:
  BardinIMDRadCorIntegrand(double z);
  ~BardinIMDRadCorIntegrand();
  double operator () (const vector<double> & x);
private:
  double fZ;
};

} // genie namespace

#endif  // _BARDIN_IMD_RADIATIVE_CORRECTIONS_PARTIAL_XSEC_H_
