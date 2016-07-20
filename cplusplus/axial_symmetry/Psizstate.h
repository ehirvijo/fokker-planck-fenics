
#ifndef __PSIZSTATE_H
#define __PSIZSTATE_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric ion psi
  // infinite space Poisson problem
  class Psizstate  : public Expression
  {
  public:
    // Evaluation routine
    void eval(Array<double>& values, const Array<double>& x) const;
    double psiz0,psi0,psif,t,t0,gamma_psi,alpha,Ti,Ti0,Tif,Tnorm,mu;
    int ionfunc;
    void compute_coeffs();

  private:
 
  };
  
}

#endif
