// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __JAC_H
#define __JAC_H

#include <dolfin.h>

namespace dolfin
{

  // A class for defining the coordinate space Jacobian
  class Jacobian : public Expression 
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
  };

}

#endif
