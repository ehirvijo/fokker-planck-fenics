
#ifndef __GDEF_H
#define __GDEF_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific type of Expression for user defined
  // source terms 
  class GSource  : public Expression
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
    double g0,gmax,t,t0,gamma_g;
    void compute_coeffs();

  private:

  };

}

#endif
