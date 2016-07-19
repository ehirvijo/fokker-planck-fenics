
#ifndef __PHIISTATE_H
#define __PHIISTATE_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric ion phi
  // infinite space Poisson problem
  class Phiistate  : public Expression
  {
  public:

    Phiistate() : t(0) {}

    // Evaluation routine
    void eval(Array<double>& values, const Array<double>& x) const;
    double phii0,phi0,phif,t,t0,gamma_phi,alpha,Ti,Ti0,Tif,Tnorm,mu;
    int ionfunc;
    void compute_coeffs();

  private:
    
  };
  
}

#endif
