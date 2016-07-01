
#include "Psiistate.h"

using namespace dolfin;

// A method to compute the value
void Psiistate::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double rx2 = x[0]*x[0]+x[1]*x[1];
 
  // Limit[Erf[x]/x,x->0]=2/sqrt(mypi) 
  if (sqrt(rx2)<1.e-6) {
    values[0] = -psii0/(8.0*pow(mypi,1.5)*alpha);
  } else {
    values[0] = -psii0*(exp(-alpha*alpha*rx2)/sqrt(mypi)+(1.0/(2*alpha*sqrt(rx2))+alpha*sqrt(rx2))*erf(sqrt(rx2)))/(8.0*mypi*alpha);
  }
  
}

// A method to compute the coefficents at this time
void Psiistate::compute_coeffs()
{
  psii0 = psi0*(1.0 + psif*(1.0 - exp(-gamma_psi*t)));
  Ti=Ti0*(1.0 + Tif*(1.0 - exp(-gamma_psi*t)));
  alpha = sqrt(Tnorm/(mu*Ti));
}

