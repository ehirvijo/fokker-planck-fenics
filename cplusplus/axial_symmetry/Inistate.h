
#ifndef __ISTATE_H
#define __ISTATE_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific type of Expression for user defined
  // initial state.
  class Inistate  : public Expression
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
  };

}

#endif
