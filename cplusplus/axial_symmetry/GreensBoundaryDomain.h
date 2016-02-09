// Copyright (C) 2016 Hirvijoki-Pfefferle
//
// This file is part of Fokker-Planck-Fenics

#ifndef __GBD_H
#define __GBD_H

#include <dolfin.h>

namespace dolfin
{

  //------------------------------------------
  // Class definition for the subdomain on which
  // we define the Green's function boundary
  // condition
  class GreensBoundaryDomain : public SubDomain 
  {
  private:

    // set the tolerance parameter used for
    // checking if the point is on the boundary
    static const double tol;
    
  public:

    // function used in fenics to flag boundary elements
    bool inside(const Array<double>& x, bool on_boundary) const;

  };  

}

#endif
