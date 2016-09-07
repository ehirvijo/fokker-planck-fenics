
#include "Vsq.h"

using namespace dolfin;

// Compute the value of the Jacobian at the requested point
// Remember that x[0]=r and x[1]=z
void Vsq::eval(Array<double>& values, const Array<double>& x) const
{ 
  values[0] = x[0]*x[0]+x[1]*x[1]; 
}
