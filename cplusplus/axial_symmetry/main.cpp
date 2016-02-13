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

  // ----------------------------------------------------------------
  // start the clock for measuring the execution time
  // ----------------------------------------------------------------

  std::clock_t start(std::clock());
  double duration;
  
  // ----------------------------------------------------------------
  // load mesh and the boundary and generate the finite element space
  // ----------------------------------------------------------------

  Mesh mesh("half_circle.xml"); 
  MeshFunction<size_t> boundaries(mesh, "half_circle_facet_region.xml");
  Forms::FunctionSpace V(mesh); 

  // --------------------------------------------------------------
  // Define the functions and derived expressions that we will need
  // --------------------------------------------------------------
  
  Function f(V); // the distribution function
  Function fprev(V); // the distribution function at previous time instance
  Function phi(V); // the first rosenbluth potential
  Function psi(V); // the second rosenbluth potential
  Source s; // the user specified source/sink
  KphiRZ Kphi; // Green's function for phi
  KpsiRZ Kpsi; // Green's function for psi
  Jacobian Jac; // Jacobian for the space 
  Inistate initial_state; // the user-specified initial state
    
  // ---------------------------------------------------------------
  // Define the parameters for the kinetic equation
  // ---------------------------------------------------------------

  Constant Efield(0.0); // the normalized electric field
  Constant nu(1.0); // the normalized collision frequency
  Constant dtau(1.0); // the normalized time step
  Constant zero(0.0); // Constant to be used as the boundary condition for the kinetic equation
  // ---------------------------------------------------------------
  // Initialize the bilinear, linear and integral forms
  // ---------------------------------------------------------------
    
  Forms::Form_aLaplace a_laplace(V, V);// bilinear volume form for the laplace operator
  Forms::Form_aProjection a_initial_projection(V, V); // bilinear volume form for the projection
  Forms::Form_aKinetic a_kinetic(V, V);// bilinear form for the time discrete kinetic equation
  Forms::Form_LProjection L_phi(V); // linear form for the Poisson equation for phi
  Forms::Form_LProjection L_psi(V); // linear form for the Poisson equation for psi
  Forms::Form_LProjection L_initial_projection(V); // linear form for the projection
  Forms::Form_LKinetic L_kinetic(V); // linear form for the time discrete kinetic equation
  Forms::Form_weightedIntegral phi_form(mesh,fprev,Kphi,Jac); // form for computing phi directly
  Forms::Form_weightedIntegral psi_form(mesh,fprev,Kpsi,Jac); // form for computing psi directly

  // ---------------------------------------------------------------
  // Assign the functions and parameters for bilinear and linear 
  // forms. These should be pass-by-reference, so updating the forms 
  // should be possible by updating the variables to which the forms 
  // are linked to. 
  // ---------------------------------------------------------------

  a_laplace.J = Jac;
  a_initial_projection.J = Jac;
  a_kinetic.J = Jac;
  a_kinetic.phi = phi;
  a_kinetic.psi = psi;
  a_kinetic.dtau = dtau;
  a_kinetic.E = Efield;
  a_kinetic.nu = nu;
  L_kinetic.J = Jac;
  L_kinetic.dtau = dtau;
  L_kinetic.S = s;
  L_kinetic.f = fprev;
  L_phi.J = Jac;
  L_psi.J = Jac;
  L_initial_projection.J = Jac;
  L_initial_projection.f = initial_state;
  L_phi.f = fprev; 
  L_psi.f = phi;
  
  // ---------------------------------------------------------------------
  // Assign the greens function solutions to the potential equations
  // ---------------------------------------------------------------------

  BCphi phi_greens_solution(&phi_form,&Kphi); // Greens function solution for phi
  BCpsi psi_greens_solution(&psi_form,&Kpsi); // Greens function solution for psi

  // ---------------------------------------------------------------------
  // Assign the boundary conditions for the equations
  // ---------------------------------------------------------------------
  DirichletBC boundary_condition_phi(V, phi_greens_solution, boundaries, 1); 
  DirichletBC boundary_condition_psi(V, psi_greens_solution, boundaries, 1); 
  DirichletBC boundary_condition_f (V, zero , boundaries, 1);

  // ---------------------------------------------------------------------
  // Project the initial state to the finite element space and take a copy
  // ---------------------------------------------------------------------
  solve( a_initial_projection == L_initial_projection, f, boundary_condition_f);
  fprev=Function(f);

  // ---------------------------------------------------------------------
  // Set up a file for storing the solution (for creating a movie with paraview)
  // ---------------------------------------------------------------------

  File file_f("f.pvd");
  File file_phi("phi.pvd");
  File file_psi("psi.pvd");


  // ---------------------------------------------------------------------
  // Loop over time, this should be extremely simple now because all 
  // functions are already assigned to the form.
  // ---------------------------------------------------------------------
  
  double t(0.0);
  int nt(40);
  for (int it = 1; it <= nt; it++) {
    
    // ---------------------------------------------------------------------
    // On the fly plotter for simple studies, use paraview to generate 
    // movies.
    // ---------------------------------------------------------------------

    std::cout<<"time step: "<<it<<"/"<<nt<<std::endl;
    //plot(f,std::string("distribution function f"));

    // ---------------------------------------------------------------------
    // Solve for the potential functions
    // ---------------------------------------------------------------------
    
    solve(a_laplace == L_phi, phi, boundary_condition_phi);
    solve(a_laplace == L_psi, psi, boundary_condition_psi);
    
    // ---------------------------------------------------------------------
    // solve for the next value of the distribution function
    // ---------------------------------------------------------------------
    
    solve(a_kinetic == L_kinetic, f, boundary_condition_f);
  
    // ---------------------------------------------------------------------
    // update the fprev from the current f
    // ---------------------------------------------------------------------

    fprev=Function(f); // Check that this doesn't cause memory leaks!!!

    // ---------------------------------------------------------------------
    // Write time stamps of f, phi, and psi to separate files.
    // ---------------------------------------------------------------------
    t= t+it*dtau;
    f.rename("f","f");
    phi.rename("phi","phi");
    psi.rename("psi","psi");
    file_f << f, t;
    file_phi << phi, t;
    file_psi << psi, t;    
	
  }
  

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout<<"duration: "<< duration <<std::endl;
  
  return 0;
}
