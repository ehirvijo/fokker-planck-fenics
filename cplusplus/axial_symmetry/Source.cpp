
#include "Source.h"

using namespace dolfin;

// Evaluate the value of the source at requested point
// Remember that the coordinates are x[0] = r, x[1] = z
void Source::eval(Array<double>& values, const Array<double>& x) const
{
  double rx2 = x[0]*x[0]+x[1]*x[1];

  values[0] = src0*exp(-rx2/srcsig);
}

// A method to compute the coefficents at this time
void Source::compute_coeffs()
{
  if (t>=t0 && t<=tf) {
    src0 = srcmax*(1.0 - exp(-gamma_src*(t-t0)));
  } else {
    src0 = 0.0;
  }
}

