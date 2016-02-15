// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "BCassembler.h"

using namespace dolfin;

// The class constructor
BCassembler::BCassembler(Forms::Form_weightedIntegral *gs,KRZ *k): _GS(gs),_K(k)
{}
  
// evaluation routine
void BCassembler::eval(Array<double>& values, const Array<double>& x) const
{
  // set up the greens function
  _K->set(x);
  _GS->k=*_K;
  values[0] = assemble(*_GS);
}
