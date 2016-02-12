#include <dolfin.h>
#include "Forms.h"
#include "BCphi.h"
#include "BCpsi.h"
#include "KphiRZ.h"
#include "KpsiRZ.h"
#include "Source.h"
#include "Inistate.h"
#include "Jacobian.h"


using namespace dolfin;

// main programm
int main()
{

  // for clocking the time
  std::clock_t start(std::clock());
  double duration;
    
  // load mesh from the file
  Mesh mesh("half_circle.xml");

  // load the boundaries, specified with "Physical Line" in the .geo file
  MeshFunction<size_t> boundaries(mesh, "half_circle_facet_region.xml");

  // Construct the Finite Element space (see ufl file)
  Forms::FunctionSpace V(mesh);

  // Define the functions we will need
  Function f(V);
  Function phi(V);
  Function psi(V);
  
  // Initialize Jacobian, bilinear and linear forms, and set the Jacobian to the forms
  Jacobian Jac;
  Forms::Form_aLaplace a_laplace(V, V);
  Forms::Form_aProjection a_initial_projection(V, V);
  Forms::Form_LProjection L_phi(V);
  Forms::Form_LProjection L_psi(V);
  Forms::Form_LProjection L_initial_projection(V);
  a_laplace.J = Jac;
  a_initial_projection.J = Jac;
  L_phi.J = Jac;
  L_psi.J = Jac;
  L_initial_projection.J=Jac;
 
  // project the initial state to the finite element basis
  Inistate f_inistate;
  DirichletBC bc_inistate(V, f_inistate, boundaries, 1); // see the .geo file 
  L_initial_projection.f = f_inistate; 
  solve( a_initial_projection == L_initial_projection, f, bc_inistate);
  
  // Solve for the potential functions
  KphiRZ Kphi; // Green's function for phi
  KpsiRZ Kpsi; // Green's function for psi
  Forms::Form_weightedIntegral phi_form(mesh,f,Kphi,Jac); 
  Forms::Form_weightedIntegral psi_form(mesh,f,Kpsi,Jac);
  BCphi phi_greens_solution(&phi_form,&Kphi); // Greens function solution for phi
  BCpsi psi_greens_solution(&psi_form,&Kpsi); // Greens function solution for psi
  L_phi.f = f; // source for the phi equation
  L_psi.f = phi; // source for the phi equation
  DirichletBC boundary_condition_phi(V, phi_greens_solution, boundaries, 1); // see the .geo file 
  DirichletBC boundary_condition_psi(V, phi_greens_solution, boundaries, 1); // see the .geo file 
  solve(a_laplace == L_phi, phi, boundary_condition_phi);
  solve(a_laplace == L_psi, psi, boundary_condition_psi);
  
  // Timing
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout<<"duration: "<< duration <<std::endl;
  
  // Plot solution
  plot(f);
  plot(phi);
  plot(psi);
  interactive();
  
  return 0;
}
