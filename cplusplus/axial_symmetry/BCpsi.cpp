// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "BCpsi.h"

using namespace dolfin;

// The class constructor
BCpsi::BCpsi(Forms::Form_I *gs,KpsiRZ *k): _GS(gs),_K(k)
{}
  
// evaluation routine
void BCpsi::eval(Array<double>& values, const Array<double>& x) const
{
  // set up the greens function
  _K->set(x);
  _GS->k=*_K;
  values[0] = assemble(*_GS);
}
