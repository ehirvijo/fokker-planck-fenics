
#ifndef __VPAR_H
#define __VPAR_H

#include <dolfin.h>

namespace dolfin
{

  // A class for defining the coordinate space Vparallel
  class Vparallel : public Expression 
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
  };

}

#endif
