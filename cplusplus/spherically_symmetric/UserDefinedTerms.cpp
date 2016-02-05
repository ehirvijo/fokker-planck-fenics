// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics
//
// Modified by David Pfefferle 2016
//
// First added:  2016-02-05

#include "UserDefinedTerms.h"

using namespace dolfin;

//-----------------------------------------------------------------------------


void Source::eval(Array<double>& values, const Array<double>& x) const
{
  values[0] = exp(-x[0]*x[0]);
}
//-----------------------------------------------------------------------------

void Jacobian::eval(Array<double>& values, const Array<double>& x) const
{ 
  values[0] = x[0]*x[0]; 
}

//----------------------------------------------------------------------

// tolerance near right boundary (should be changed)
const double VeloBoundary::tol=1.e-15;

// function identifying the boundary
bool VeloBoundary::inside(const Array<double>& x, bool on_boundary) const
{  return on_boundary and abs(x[0]-rightbound)<tol;}
  
// constructor given the upper limit (again, not generic)
VeloBoundary::VeloBoundary(const double& ul) : rightbound(ul){}

//-----------------------------------------------------------------------------

// The logarithmically space mesh constructor 
LogarithmicIntervalMesh::LogarithmicIntervalMesh(const size_t& dim,const double& leftb, const double& rightb) : IntervalMesh(dim,leftb,rightb)
{
  double loglow(-2.);   // you want to make this more general
  // increment in log scale
  double pivot((log10(rightb)-loglow)/double(dim-1));
  //define mesh elements from 1 to dim-1
  for(size_t i=1;i<dim;i++){
    this->coordinates()[i] = pow(10,loglow + double(i-1)*pivot); // check this
  }
}
