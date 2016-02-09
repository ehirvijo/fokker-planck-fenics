#include <dolfin.h>
#include "Forms.h"
#include "GreensBoundaryCondition.h"
#include "GreensBoundaryDomain.h"
#include "KphiRZ.h"
#include "Source.h"
#include "Jacobian.h"


using namespace dolfin;

// main programm
int main(int argc,char** argv)
{
  std::clock_t start(std::clock());
  double duration;
  
  // parameters for mesh creation (could appear as argv, argc)
  size_t nr = 51;
  size_t nz = 50;
  double rmin = 0. ;
  double rmax = 10. ;
  double zabsmax = 10. ;
  
  if(argc==6) {
    nr=std::atoi(argv[1]);
    nz=std::atoi(argv[2]);
    rmin=std::atof(argv[3]);
    rmax=std::atof(argv[4]);
    zabsmax=std::atof(argv[5]);

    std::cout<<"using nr="<<nr<<", nz="<<nz<<", rmin="<<rmin<<", rmax="<<rmax<<", zabsmax="<<zabsmax<<std::endl;
  }
  
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

  // Timing
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout<<"duration: "<< duration <<std::endl;
  
  // Plot solution
  plot(u);
  interactive();

  return 0;
}
