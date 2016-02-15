// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __KRZ_H
#define __KRZ_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric Greens function for the 
  // infinite space Poisson problem
  class KRZ  : public Expression
  {
    
    // protected (so that they are passed on to children class, but
    // not publically accessible)
  protected:
    double r0;
    double z0;
    
  public:
    
    //virtual void eval(Array<double>& values, const Array<double>& x)
    
    // const; routine to set up the parameters
    void set(const Array<double>& x);
    
  };
  
}

#endif
