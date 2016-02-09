// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#include "GreensBoundaryDomain.h"

using namespace dolfin;

// Indentify the boundary domain for the Dirichlet condition
// Remember that the x[0]=r and x[1]=z. The boundary domain
// for the Dirichlet condition is the part of boundary where
// r>0, i.e. the boundary domain r==0 is not included.
//
bool GreensBoundaryDomain::inside(const Array<double>& x, bool on_boundary) const
{  
  return on_boundary and abs(x[0])>tol;
}

const double GreensBoundaryDomain::tol=1.e-15;
