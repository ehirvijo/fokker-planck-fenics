
#include "Psizstate.h"

using namespace dolfin;

// A method to compute the value
void Psizstate::eval(Array<double>& values, const Array<double>& x) const
{
  double mypi(3.14159265359);
  double rx2 = x[0]*x[0]+x[1]*x[1];
 
  // Limit[Erf[x]/x,x->0]=2/sqrt(mypi) 
  if (sqrt(rx2)<1.e-6) {
    values[0] = -psiz0/(8.0*pow(mypi,1.5)*alpha);
  } else {
    values[0] = -psiz0*(exp(-alpha*alpha*rx2)/sqrt(mypi)+(1.0/(2*alpha*sqrt(rx2))+alpha*sqrt(rx2))*erf(alpha*sqrt(rx2)))/(8.0*mypi*alpha);
  }
  
}

// A method to compute the coefficents at this time
void Psizstate::compute_coeffs()
{
  if (ionfunc==0) {
    if (t>=t0 && t0!=-1) {
      psiz0 = psi0*(1.0 + psif*(1.0 - exp(-gamma_psi*(t-t0))));
      Ti = Ti0*(1.0 + Tif*(1.0 - exp(-gamma_psi*(t-t0))));
    } else {
      psiz0 = psi0;
      Ti = Ti0;
    }
  } else {
    psiz0=psi0;
    Ti=Ti0;
  } 
  alpha = sqrt(Tnorm/(mu*Ti));
}

