
#ifndef __VSQ_H
#define __VSQ_H

#include <dolfin.h>

namespace dolfin
{

  // A class for defining the coordinate space Vparallel
  class Vsq : public Expression 
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
  };

}

#endif
