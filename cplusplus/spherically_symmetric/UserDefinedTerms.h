// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics
//
// Modified by David Pfefferle 2016
//
// First added:  2016-02-05

#ifndef __UDEF_H
#define __UDEF_H

#include <dolfin.h>
#include "SphericallySymmetric.h"

namespace dolfin
{
  
  class VeloBoundary : public SubDomain
  {
  private:
    static const double tol;
    double rightbound;              // memory of right boundary point, this is not generic, should find a better way
    
  public:
    // function used in fenics to flag boundary elements
    bool inside(const Array<double>& x, bool on_boundary) const;
    
    // constructor given the upper limit (again, not generic, please think of something better)
    explicit VeloBoundary(const double& ul);
  };  
  
  // Source term (right-hand side)
  class Source : public Expression
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
  };
  
  // Jacobian term
  class Jacobian : public Expression
  {
  public:
    void eval(Array<double>& values, const Array<double>& x) const;
  };
  
  // The logarithmically spaced mesh
  class LogarithmicIntervalMesh : public IntervalMesh
  {
  public:
    // explicit constructor
    explicit LogarithmicIntervalMesh(const size_t&,const double&,const double&);
  };
  
  class GreenBC : public Expression
  {
  public:
    SphericallySymmetric::Form_I B;
    
    void eval(Array<double>& values, const Array<double>& x) const;

explicit GreenBC(const Mesh& mesh, const Expression& f,const Expression& jac);
  };

}

#endif
