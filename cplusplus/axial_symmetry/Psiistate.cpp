
#include "Psiistate.h"

using namespace dolfin;

// A method to compute the value
void Psiistate::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double rx2 = x[0]*x[0]+x[1]*x[1];
 
  // Limit[Erf[x]/x,x->0]=2/sqrt(mypi) 
  if (sqrt(rx2)<1.e-6) {
    values[0] = -psii0/4.;
  } else {
    values[0] = -psii0*(2*sqrt(rx2)*exp(-rx2)+sqrt(mypi)*(1+2*rx2)*erf(sqrt(rx2)))/(16*sqrt(rx2));
  }
  
}

// A method to compute the coefficents at this time
void Psiistate::compute_coeffs()
{
  psii0 = psimax*((1.0 + exp(-gamma_psi*t))/2.0);
}

