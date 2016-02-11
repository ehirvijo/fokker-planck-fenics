#include <dolfin.h>
#include "Forms.h"
#include "BCphi.h"
#include "BCpsi.h"
#include "GreensBoundaryDomain.h"
#include "KphiRZ.h"
#include "KpsiRZ.h"
#include "Source.h"
#include "Jacobian.h"


using namespace dolfin;

// main programm
int main()
{
  std::clock_t start(std::clock());
  double duration;
    
  // load mesh from the file
  Mesh mesh("half_circle.xml");

  // Define Finite Element space (see ufl file)
  Forms::FunctionSpace V(mesh);
  
  // define the boundary on which
  // the Dirichlet condition will
  // beapplied
  GreensBoundaryDomain boundary;

  // Define the expressions that are common for
  // both potential equations
  Jacobian J;
  Source f;
  Forms::BilinearForm a(V, V);
  Forms::LinearForm L(V);
  L.J = J;
  a.J = J;

  //------------------------------------------------------
  // Set up the stuff for the potential phi
  //
  KphiRZ Kphi;
  L.f = f;
  Forms::Form_I gsphi(mesh,f,Kphi,J);
  BCphi gbcphi(&gsphi,&Kphi);
  DirichletBC bcphi(V, gbcphi, boundary);
  
  // Compute solution
  Function uphi(V);
  solve(a == L, uphi, bcphi);
  
  //------------------------------------------------------
  // set up the stuff for the second potential
  KpsiRZ Kpsi;
  L.f = uphi;
  Forms::Form_I gspsi(mesh,f,Kpsi,J);
  BCpsi gbcpsi(&gspsi,&Kpsi);
  DirichletBC bcpsi(V, gbcpsi, boundary);
  
  // Compute solution
  Function upsi(V);
  solve(a == L, upsi, bcpsi);

  // // Save solution in VTK format
  // File file("axially_symmetric.pvd");
  // file << u;

  // Timing
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout<<"duration: "<< duration <<std::endl;
  
  // Plot solution
  plot(uphi);
  plot(upsi);
  interactive();

  return 0;
}
