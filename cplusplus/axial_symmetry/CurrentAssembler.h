// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __CURA_H
#define __CURA_H

#include <dolfin.h>
#include "Forms.h"

namespace dolfin
{
  
  //--------------------------------------------------
  // Class to define the expression for the Dirichlet
  // condition at the "GreensBoundaryDomain"
  class CurrentAssembler : public Expression
  {
  public:
    
    // Stores a reference to a UFL form
    Forms::Form_LProjection *_IG;
    
    // a routine to evaluate the current moment
    void eval(Array<double>& values, const Array<double>& x) const;
    
    // Constructor, assigns the Expressions "f"
    // (the source in the Poisson equation), 
    // "k" (the Green's function), and "j" 
    // (the Jacobian of the coordinate system) 
    // to the form "B", which then can be evaluated
    // by calling the "eval" function
    explicit BCassembler(Forms::Form_weightedIntegral *gs, KRZ *k);
    
  };
  
}

#endif
