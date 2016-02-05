#include <dolfin.h>
#include "UserDefinedTerms.h"
#include "SphericallySymmetric.h"

using namespace dolfin;

// main programm
int main()
{
  // parameters for mesh creation (could appear as argv, argc)
  size_t nx=30;
  double xmin=0.01;
  double xmax=10;
  
  // Create mesh
  LogarithmicIntervalMesh mesh(nx-1,xmin,xmax);
  
  // Define Finite Element space (see ufl file)
  SphericallySymmetric::FunctionSpace V(mesh);
  
  // boundary object
  VeloBoundary boundary(xmax);
  // boundary condition
  //constant is evaluated as TruncGreenIntegral(xmin,xmax)
  Constant u0(0.0); // zero for the moment

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
  File file("spherically_symmetric.pvd");
  file << u;

  // Plot solution
  plot(u);
  interactive();

  return 0;
}
