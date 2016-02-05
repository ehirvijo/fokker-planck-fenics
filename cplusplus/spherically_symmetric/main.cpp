#include <dolfin.h>
#include "UserDefinedTerms.h"
#include "SphericallySymmetric.h"

using namespace dolfin;

// main programm
int main()
{
  // parameters for mesh creation (could appear as argv, argc)
  size_t nx=30;
  double xmin=0;
  double xmax=10;
  
  // Create mesh
  LogarithmicIntervalMesh mesh(nx-1,xmin,xmax);
  
  // Define function space for the mesh using Continuous Galerkin
  // (Lagrange) functions of order p on each element
  SphericallySymmetric::FunctionSpace V(mesh);
  
  // boundary object
  VeloBoundary boundary(xmax);
  // boundary condition
  //constant is evaluated as TruncGreenIntegral(xmin,xmax)
  Constant u0(0.0);

  DirichletBC bc(V, u0, boundary);

  Source f;
  Jacobian x2;

  // Define variational problem
  SphericallySymmetric::BilinearForm a(V, V);
  SphericallySymmetric::LinearForm L(V);
  L.f = f;
  L.J = x2;
  a.J = x2;

  // Compute solution
  Function u(V);
  solve(a == L, u, bc);

  // Save solution in VTK format
  File file("biharmonic.pvd");
  file << u;

  // Plot solution
  plot(u);
  interactive();

  return 0;
}
