// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __KPSIRZ_H
#define __KPSIRZ_H

#include <dolfin.h>
#include "KRZ.h"

namespace dolfin
{
  
  //-----------------------------------------------
  // A specific expression for defining the 
  // axially symmetric Greens function for the 
  // infinite space biharmonic problem
  class KpsiRZ  : public KRZ
  {
    // Evaluation routine
    void eval(Array<double>& values, const Array<double>& x) const;
    
  };
  
}

#endif
