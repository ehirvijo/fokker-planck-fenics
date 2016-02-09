// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "GreensBoundaryCondition.h"

using namespace dolfin;

// The class constructor
GreensBoundaryCondition::GreensBoundaryCondition(const Mesh& m, const Expression& f, const Expression& k, const Expression& j) : greensSolution(Forms::Form_I(m,f,k,j))
{}
  
// evaluation routine
void GreensBoundaryCondition::eval(Array<double>& values, const Array<double>& x) const
{
  // set up the greens function
  K.set(x);
  greensSolution.k=K;
  values[0] = assemble(greensSolution);
}
