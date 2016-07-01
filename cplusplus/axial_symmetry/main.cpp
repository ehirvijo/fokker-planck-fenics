// 2016 Hirvijoki-Pfefferle-Brennan Princeton Plasma Physics Lab
#include <dolfin.h>
#include <iostream>
#include <fstream>
#include <string>
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
#include "Vparallel.h"

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
  Vparallel Vpar; // "function" for assembling the plasma current form
  Jacobian Jac; // Jacobian for the space 
  Psiistate psii_state; // the user-specified ion psi state
  Phiistate phii_state; // the user-specified ion phi state
  Inistate initial_state; // the user-specified initial state
    
  // ---------------------------------------------------------------
  // Define the parameters for the kinetic equation
  // ---------------------------------------------------------------

  Constant Efield(0.005); // the normalized electric field
  Constant nu(4.0e-6); // the normalized collision frequency
  Constant nui(1.6e-5); // the normalized collision frequency ions = nu*Zeff^2
  Constant dtau(1.0); // the normalized time step in slowing times
  Constant mu(2.72443712e-04); // the electron to ion mass ratio 9.1093898e-31/3.3435860e-27
  Constant zero(0.0); // Constant to be used as the boundary condition for the kinetic equation
  double Tnorm(1.0e4); // Normalization temperature for vo
  double Ti0(1.0e4); // Initial ion temperature
  double Tif(-0.9); // Fraction of ion temperature change in collapse time
  double ni0(1.0); // Initial ion density fraction to reference density
  double nif(0.5); // fraction of ion density change in collapse time
  double gamma_c(0.3); // The thermal collapse rate (inverse timescale)
  double gamma_g(0.3); // The solution dependent source rise rate 
  double gamma_src(0.3); // The source rise rate
  double gmax(0.0); // The solution dependent source factor
  double srcmax(0.0); // The source max
  int nt(1000);
  
  // ---------------------------------------------------------------
  // Read the input.dat file.
  // Values are assumed specified as "variable_name value" on each line
  // ---------------------------------------------------------------
  std::string line;
  std::ifstream file_inp("input.dat");
  float val;
  if (file_inp.is_open()) {
    while(getline(file_inp,line)) {
      std::istringstream in(line);
      std::string(type);
      in>>type;
      if(type == "srcmax") {
        in>>srcmax;
      } else if (type == "gamma_src") {
        in>>gamma_src;
      } else if (type == "nt") {
        in>>nt;
      } else if (type == "Tnorm") {
        in>>Tnorm;
      } else if (type == "Ti0") {
        in>>Ti0;
      } else if (type == "Tif") {
        in>>Tif;
      } else if (type == "Efield") {
        in>>val;
        Efield=val;
      } else if (type == "dtau") {
        in>>val;
        dtau=val;
      } else if (type == "nu") {
        in>>val;
        nu=val;
      } else if (type == "nui") {
        in>>val;
        nui=val;
      } else {
        std::cout << "Unknown setting input " <<type<<std::endl;
      }
    }
    file_inp.close();
  } else {
    std::cout << "No input file" << std::endl;
  }

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
  Forms::Form_weightedIntegral current_form(mesh,fprev,Vpar,Jac); // form for computing psi directly

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
  
  current_form.f=fprev;
  current_form.k=Vpar;
  current_form.J=Jac;

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
  // Set ion potential parameters
  // the redundant variables shoudl be collapsed into one
  // ---------------------------------------------------------------------
  phii_state.phi0=ni0;
  psii_state.psi0=ni0;
  phii_state.phif=nif;
  psii_state.psif=nif;
  phii_state.Ti0=Ti0; 
  psii_state.Ti0=Ti0; 
  phii_state.Tif=Tif;
  psii_state.Tif=Tif;
  phii_state.Tnorm=Tnorm;
  psii_state.Tnorm=Tnorm;
  phii_state.mu=mu;
  psii_state.mu=mu;
  phii_state.gamma_phi=gamma_c;
  psii_state.gamma_psi=gamma_c;
  phii_state.compute_coeffs();
  psii_state.compute_coeffs();
  phii.interpolate(phii_state); 
  psii.interpolate(psii_state); 

  // ---------------------------------------------------------------------
  // Set up solution independent source
  // ---------------------------------------------------------------------
  s.gamma_src=gamma_src;
  s.srcmax=srcmax;
  s.compute_coeffs();

  // ---------------------------------------------------------------------
  // Set up solution dependent source
  // ---------------------------------------------------------------------
  gs.gamma_g=gamma_g;
  gs.gmax=gmax;
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
  std::ofstream file_dis("discharge.dat");
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
  double jpar(0.0);
  for (int it = 1; it <= nt; it++) {

    std::cout<<"time step: "<<it<<"/"<<nt<<" time: "<<t<< " jpar: "<<jpar<<std::endl;

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
    // Calculate the parallel current density for the electric field update
    // ---------------------------------------------------------------------
    jpar=assemble(current_form);

    // ---------------------------------------------------------------------
    // Advance the time
    // ---------------------------------------------------------------------
    t= t+dtau;

    // ---------------------------------------------------------------------
    // Write time stamps of f, phi, psi etc. to separate files. 
    // ---------------------------------------------------------------------
    if (it % 10==0) {
      file_f << f, t;
      file_phi << phi, t;
      file_psi << psi, t;    
      file_phii << phii, t;
      file_psii << psii, t;    
    };

    file_dis << it <<" "<< t <<" "<< dtau <<" "<< jpar <<" "<< Efield 
             << " "<< phii_state.Ti <<" "<< phii_state.phii0 << "\n";
	
  }
  swatch.stop();

  list_timings();
  
  file_dis.close();

  interactive();
  return 0;
}


