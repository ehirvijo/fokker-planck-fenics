// 2016 Hirvijoki-Pfefferle-Brennan Princeton Plasma Physics Lab
#include <dolfin.h>
#include "Forms.h"
#include "BCassembler.h"
#include "KphiRZ.h"
#include "KpsiRZ.h"
#include "Source.h"
#include "GSource.h"
#include "Inistate.h"
#include "Phiistate.h"
#include "Psiistate.h"
#include "Jacobian.h"

using namespace dolfin;

// main programm
int main()
{

  // ----------------------------------------------------------------
  // start the clock for measuring the execution time
  // ----------------------------------------------------------------

  Timer swatch("Initialisation, mesh loading and assigning");
  
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
  Function phii(V); // the first rosenbluth potential ions
  Function psi(V); // the second rosenbluth potential
  Function psii(V); // the second rosenbluth potential ions
  Function g(V); // the parallel dependent sink function on the grid
  GSource gs; // the parallel dependent expression
  Source s; // the user specified source/sink expression
  KphiRZ Kphi; // Green's function for phi
  KpsiRZ Kpsi; // Green's function for psi
  Jacobian Jac; // Jacobian for the space 
  Psiistate psii_state; // the user-specified ion psi state
  Phiistate phii_state; // the user-specified ion phi state
  Inistate initial_state; // the user-specified initial state
    
  // ---------------------------------------------------------------
  // Define the parameters for the kinetic equation
  // ---------------------------------------------------------------

  Constant Efield(0.0); // the normalized electric field
  Constant nu(1.0); // the normalized collision frequency
  Constant nui(1.0); // the normalized collision frequency ions = nu*Zeff^2
  Constant dtau(1.0); // the normalized time step
  Constant mu(2.72443712e-04); // the electron to ion mass ratio 9.1093898e-31/3.3435860e-27
  Constant zero(0.0); // Constant to be used as the boundary condition for the kinetic equation

  // ---------------------------------------------------------------
  // Initialize the bilinear, linear and integral forms
  // ---------------------------------------------------------------
    
  Forms::Form_aLaplace a_laplace(V, V);// bilinear volume form for the laplace operator
  Forms::Form_aKinetic a_kinetic(V, V);// bilinear form for the time discrete kinetic equation
  Forms::Form_LProjection L_phi(V); // linear form for the Poisson equation for phi
  Forms::Form_LProjection L_psi(V); // linear form for the Poisson equation for psi
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
  a_kinetic.J = Jac;
  
  a_kinetic.phi = phi;
  a_kinetic.psi = psi;
  a_kinetic.phii = phii;
  a_kinetic.psii = psii;
  a_kinetic.dtau = dtau;
  a_kinetic.E = Efield;
  a_kinetic.nu = nu;
  a_kinetic.nui = nui;
  a_kinetic.mu = mu;
  a_kinetic.g = g;
  
  L_kinetic.J = Jac;
  L_kinetic.dtau = dtau;
  L_kinetic.S = s;
  L_kinetic.f = fprev;
  
  L_phi.J = Jac; 
  L_psi.J = Jac;
  L_phi.f = fprev; 
  L_psi.f = phi;
  
  // ---------------------------------------------------------------------
  // Assign the greens function solutions to the potential equations
  // ---------------------------------------------------------------------

  BCassembler phi_greens_solution(&phi_form,&Kphi); // Greens function solution for phi
  BCassembler psi_greens_solution(&psi_form,&Kpsi); // Greens function solution for psi

  // ---------------------------------------------------------------------
  // Assign the boundary conditions for the equations
  // ---------------------------------------------------------------------
  DirichletBC boundary_condition_phi(V, phi_greens_solution, boundaries, 1); 
  DirichletBC boundary_condition_psi(V, psi_greens_solution, boundaries, 1); 
  DirichletBC boundary_condition_f (V, zero , boundaries, 1);

  // ---------------------------------------------------------------------
  // Set initial ion potentials
  // ---------------------------------------------------------------------
  phii_state.phimax=1.0;
  psii_state.psimax=1.0;
  phii_state.gamma_phi=0.3;
  psii_state.gamma_psi=0.3;
  phii_state.compute_coeffs();
  psii_state.compute_coeffs();
  phii.interpolate(phii_state); 
  psii.interpolate(psii_state); 

  // ---------------------------------------------------------------------
  // Set up solution independent source
  // ---------------------------------------------------------------------
  s.gamma_src=0.3;
  s.srcmax=0.0;
  s.compute_coeffs();

  // ---------------------------------------------------------------------
  // Set up solution dependent source
  // ---------------------------------------------------------------------
  gs.gamma_g=0.3;
  gs.gmax=0.1;
  gs.compute_coeffs();
  g.interpolate(gs);

  // ---------------------------------------------------------------------
  // Project the initial state to the finite element space and take a copy
  // ---------------------------------------------------------------------
  f.interpolate(initial_state);
  fprev=f;

  // ---------------------------------------------------------------------
  // Set up a file for storing the solution (for creating a movie with paraview)
  // The rename calls ensure the variables are named correctly in the files.
  // ---------------------------------------------------------------------

  File file_f("f.pvd");
  File file_phi("phi.pvd");
  File file_phii("phii.pvd");
  File file_psi("psi.pvd");
  File file_psii("psii.pvd");
  f.rename("f","f");
  phi.rename("phi","phi"); 
  psi.rename("psi","psi");
  phii.rename("phi","phi");
  psii.rename("psi","psi");

  // ---------------------------------------------------------------------
  // Loop over time, this should be extremely simple now because all 
  // functions are already assigned to the form.
  // ---------------------------------------------------------------------

  swatch.stop();
  
  swatch = Timer("Time loop and problem solving");
  
  double t(0.0);
  int nt(40);
  for (int it = 1; it <= nt; it++) {

    std::cout<<"time step: "<<it<<"/"<<nt<<" time: "<<t<<std::endl;

    // ---------------------------------------------------------------------
    // On the fly plotting for simple studies, use paraview to generate 
    // movies.
    // ---------------------------------------------------------------------
    plot(f,std::string("distribution function f"));
    //plot(g,std::string("sink function operator"));
    //plot(phii,std::string("phii_state"));
    //plot(psii,std::string("psii_state"));
    //interactive();

    // ---------------------------------------------------------------------
    // Set ion potentials
    // ---------------------------------------------------------------------
    phii_state.t=t;
    psii_state.t=t;
    phii_state.compute_coeffs();
    psii_state.compute_coeffs();
    phii.interpolate(phii_state); 
    psii.interpolate(psii_state); 
  
    // ---------------------------------------------------------------------
    // Set the solution independent source
    // ---------------------------------------------------------------------
    s.t=t;
    s.compute_coeffs();
  
    // ---------------------------------------------------------------------
    // Set the solution dependent source
    // ---------------------------------------------------------------------
    gs.t=t;
    gs.compute_coeffs();
    g.interpolate(gs);

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

    fprev=f; // Check that this doesn't cause memory leaks!!!
    
    // ---------------------------------------------------------------------
    // Advance the time
    // ---------------------------------------------------------------------

    t= t+dtau;

    // ---------------------------------------------------------------------
    // Write time stamps of f, phi, psi etc. to separate files.  The 
    // rename calls ensure the variables are named correctly in the files.
    // ---------------------------------------------------------------------

    file_f << f, t;
    file_phi << phi, t;
    file_psi << psi, t;    
    file_phii << phii, t;
    file_psii << psii, t;    
	
  }
  swatch.stop();

  list_timings();
  
  return 0;
}


