// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "KRZ.h"

using namespace dolfin;

// A method for setting the parameters
void KRZ::set(const Array<double>& x)
{
  r0=x[0];
  z0=x[1];
}
