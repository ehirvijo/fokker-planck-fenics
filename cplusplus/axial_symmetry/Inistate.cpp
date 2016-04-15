
#include "Inistate.h"

using namespace dolfin;

// Evaluate the value of the source at requested point
// Remember that the coordinates are x[0] = r, x[1] = z
void Inistate::eval(Array<double>& values, const Array<double>& x) const
{
  // values[0] = exp(-x[0]*x[0]-(x[1]-2.0)*(x[1]-2.0)) + exp(-x[0]*x[0]-(x[1]+2.0)*(x[1]+2.0));
  values[0] = exp(-x[0]*x[0]-x[1]*x[1]);
}
