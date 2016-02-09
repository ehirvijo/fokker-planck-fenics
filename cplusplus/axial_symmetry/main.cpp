#include <dolfin.h>
#include "Forms.h"
#include "GreensBoundaryCondition.h"
#include "GreensBoundaryDomain.h"
#include "KphiRZ.h"
#include "Source.h"
#include "Jacobian.h"


using namespace dolfin;

// main programm
int main()
{
  // parameters for mesh creation (could appear as argv, argc)
  size_t nr = 31 ; 
  size_t nz = 21 ;
  double rmin = 0. ;
  double rmax = 5. ;
  double zabsmax = -5. ;
  
  // Create a simple rectangular mesh
  RectangleMesh mesh(rmin,-zabsmax,rmax,zabsmax,nr,nz, std::string("right"));

  // plot the mesh
  // plot(mesh);
  // interactive();

  // Define Finite Element space (see ufl file)
  Forms::FunctionSpace V(mesh);
  
  // define the boundary on which
  // the Dirichlet condition will
  // beapplied
  GreensBoundaryDomain boundary;

  // Define the expressions that are needed
  // to set up the Greens function boundary 
  // condition at the Dirichlet boundary
  Source f;
  Jacobian J;
  KphiRZ Kphi;
  Forms::Form_I gs(mesh,f,Kphi,J);

  // Define the variational problem
  Forms::BilinearForm a(V, V);
  Forms::LinearForm L(V);
  L.f = f;
  L.J = J;
  a.J = J;
  
  // set up the Greens function boundary condition
  // using the expressions f, K_phi, J and the mesh
  GreensBoundaryCondition gbc(&gs,&Kphi);
  DirichletBC bc(V, gbc, boundary);
  
  // Compute solution
  Function u(V);
  solve(a == L, u, bc);
  
  // Save solution in VTK format
  File file("axially_symmetric.pvd");
  file << u;

  // Plot solution
  plot(u);
  interactive();

  return 0;
}
