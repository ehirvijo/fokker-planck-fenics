#include <dolfin.h>
#include "Forms.h"
#include "BCphi.h"
#include "BCpsi.h"
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

  // load the boundaries, specified with "Physical Line" in the .geo file
  MeshFunction<size_t> boundaries(mesh, "half_circle_facet_region.xml");

  // Construct the Finite Element space (see ufl file)
  Forms::FunctionSpace V(mesh);
  
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
  DirichletBC bcphi(V, gbcphi, boundaries,1); // 1 refers to the .geo file 

  // Compute solution
  Function uphi(V);
  solve(a == L, uphi, bcphi);
  
  //------------------------------------------------------
  // set up the stuff for the second potential
  KpsiRZ Kpsi;
  L.f = uphi;
  Forms::Form_I gspsi(mesh,f,Kpsi,J);
  BCpsi gbcpsi(&gspsi,&Kpsi);
  DirichletBC bcpsi(V, gbcpsi, boundaries,1); // 1 refers to the .geo file 
    
  // Compute solution
  Function upsi(V);
  solve(a == L, upsi, bcpsi);

  // Save solution in VTK format
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
