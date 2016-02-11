// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __BCPSI_H
#define __BCPSI_H

#include <dolfin.h>
#include "KpsiRZ.h"
#include "Forms.h"

namespace dolfin
{
  
  //--------------------------------------------------
  // Class to define the expression for the Dirichlet
  // condition at the "GreensBoundaryDomain"
  class BCpsi : public Expression
  {
  public:
    
    // Stores a reference to a UFL form
    Forms::Form_I *_GS;
    KpsiRZ *_K;
    
    // a routine to evaluate the value of the Greens's function
    // solution at point x
    void eval(Array<double>& values, const Array<double>& x) const;
    
    // Constructor, assigns the Expressions "f"
    // (the source in the Poisson equation), 
    // "k" (the Green's function), and "j" 
    // (the Jacobian of the coordinate system) 
    // to the form "B", which then can be evaluated
    // by calling the "eval" function
    explicit BCpsi(Forms::Form_I *gs, KpsiRZ *k);

  };

}

#endif
