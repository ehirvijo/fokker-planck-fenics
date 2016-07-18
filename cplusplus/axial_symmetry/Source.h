
#ifndef __UDEF_H
#define __UDEF_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific type of Expression for user defined
  // source terms 
  class Source  : public Expression
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
    double src0,srcmax,srcsig,t,t0,tf,gamma_src;
    void compute_coeffs();

  private:

  };

}

#endif
