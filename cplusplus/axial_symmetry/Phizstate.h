
#ifndef __PHIZSTATE_H
#define __PHIZSTATE_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric ion phi
  // infinite space Poisson problem
  class Phizstate  : public Expression
  {
  public:

    Phizstate() : t(0) {}

    // Evaluation routine
    void eval(Array<double>& values, const Array<double>& x) const;
    double phiz0,phi0,phif,t,t0,gamma_phi,alpha,Ti,Ti0,Tif,Tnorm,mu;
    int ionfunc;
    void compute_coeffs();

  private:
    
  };
  
}

#endif
