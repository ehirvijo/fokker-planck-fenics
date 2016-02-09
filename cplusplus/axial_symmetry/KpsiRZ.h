// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __KPSIRZ_H
#define __KPSIRZ_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric Greens function for the 
  // infinite space biharmonic problem
  class KpsiRZ  : public Expression
  {
    
  private:
    
    double r0;
    double z0;
    
  public:
    
    // Evaluation routine
    void eval(Array<double>& values, const Array<double>& x) const;
    
    // routine to set up the parameters
    void set(const Array<double>& x);
    
  };
  
}

#endif
