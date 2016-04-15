
#include "CurrentAssembler.h"

using namespace dolfin;

// The class constructor
CurrentAssembler::CurrentAssembler(Forms::Form_LProjection *ig): _IG(ig)
{}
  
// evaluation routine
void CurrentAssembler::eval(Array<double>& values, const Array<double>& x) const
{
  _IG->k=*_K;
  values[0] = assemble(*_IG);
}
