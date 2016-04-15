// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __PSIISTATE_H
#define __PSIISTATE_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric ion psi
  // infinite space Poisson problem
  class Psiistate  : public Expression
  {
  public:
    // Evaluation routine
    void eval(Array<double>& values, const Array<double>& x) const;
    double psii0,psimax,t,gamma_psi;
    void compute_coeffs();

  private:
 
  };
  
}

#endif
