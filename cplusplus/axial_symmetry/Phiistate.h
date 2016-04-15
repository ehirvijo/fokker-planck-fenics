// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __PHIISTATE_H
#define __PHIISTATE_H

#include <dolfin.h>

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric ion phi
  // infinite space Poisson problem
  class Phiistate  : public Expression
  {
  public:

    Phiistate() : t(0) {}

    // Evaluation routine
    void eval(Array<double>& values, const Array<double>& x) const;
    double phii0,phimax,t,gamma_phi;
    void compute_coeffs();

  private:
    
  };
  
}

#endif
