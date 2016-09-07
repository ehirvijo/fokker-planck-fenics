
#include "Phiistate.h"

using namespace dolfin;

// A method to compute the value
void Phiistate::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double rx2 = x[0]*x[0]+x[1]*x[1];

  // Limit[Erf[x]/x,x->0]=2/sqrt(mypi)
  if (sqrt(rx2)<1.e-6) { 
    values[0] = -phii0/(alpha*sqrt(mypi));
  } else {
    values[0] = -phii0*erf(alpha*sqrt(rx2))/(2.0*pow(alpha,2)*sqrt(rx2));
  }
  
}

// A method to compute the coefficents at this time
void Phiistate::compute_coeffs() 
{
  if (ionfunc==0) {
    if (t>=t0) {
      phii0 = phi0*(1.0 + phif*(1.0 - exp(-gamma_phi*(t-t0))));
      Ti = Ti0*(1.0 + Tif*(1.0 - exp(-gamma_phi*(t-t0))));
    } else {
      phii0 = phi0;
      Ti = Ti0;
    }
  } else {
    phii0 = phi0;
    Ti = Ti0;
  }
  alpha = sqrt(Tnorm/(mu*Ti));
}
