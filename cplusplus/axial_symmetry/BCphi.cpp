// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "BCphi.h"

using namespace dolfin;

// The class constructor
BCphi::BCphi(Forms::Form_weightedIntegral *gs,KphiRZ *k): _GS(gs),_K(k)
{}
  
// evaluation routine
void BCphi::eval(Array<double>& values, const Array<double>& x) const
{
  // set up the greens function
  _K->set(x);
  _GS->k=*_K;
  values[0] = assemble(*_GS);
}
