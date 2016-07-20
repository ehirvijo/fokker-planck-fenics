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
#include "Phizstate.h"
#include "Psizstate.h"
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
  
  Function f(V); // distribution function
  Function fprev(V); // distribution function at previous time instance
  Function phi(V); // first rosenbluth potential
  Function phii(V); // first rosenbluth potential ions
  Function phiz(V); // first rosenbluth potential impurity ions
  Function psi(V); // second rosenbluth potential
  Function psii(V); // second rosenbluth potential ions
  Function psiz(V); // second rosenbluth potential impurity ions
  Function g(V); // parallel dependent sink function on the grid
  GSource gs; // parallel dependent expression
  Source s; // user specified source/sink expression
  KphiRZ Kphi; // Green's function for phi
  KpsiRZ Kpsi; // Green's function for psi
  Vparallel Vpar; // "function" for assembling the plasma current form
  Jacobian Jac; // Jacobian for the space 
  Psiistate psii_state; // user-specified ion psi state
  Phiistate phii_state; // user-specified ion phi state
  Psizstate psiz_state; // user-specified impurity ion psi state
  Phizstate phiz_state; // user-specified impurity ion phi state
  Inistate initial_state; // user-specified initial state
    
  // ---------------------------------------------------------------
  // Define the parameters for the kinetic equation
  // ---------------------------------------------------------------

  Constant Efield(0.005); // normalized electric field
  Constant nu(1.0); // normalized collision frequency
  Constant nui(1.0); // normalized collision frequency ions = nu*Z0
  Constant nuz(18.0); // normalized collision frequency impurity ions = nu*Zi
  Constant dtau(1.0); // normalized time step in slowing times
  Constant mu(2.72443712e-04); // electron to ion mass ratio 9.1093898e-31/3.3435860e-27
  Constant muz(6.81995875e-06); // electron to impurity ion mass ratio 9.1093898e-31/3.3435860e-27/39.948  Argon assumed for Zi and mass here
  Constant zero(0.0); // constant to be used as the boundary condition for the kinetic equation
  Constant one(1.0); // constant to be used integrate the density
  double Tnorm(1.0e4); // normalization temperature for vo

  // ---------------------------------------------------------------
  // Define the parameters controlling the run scenario
  // ---------------------------------------------------------------

  double Ti0(1.0e4); // initial ion temperature
  double Tiz(1.0e4); // initial impurity ion temperature
  double Tif(0.0); // fraction of ion temperature change in collapse time
  double ni0(1.0); // initial ion density fraction to reference density
  double nif(0.0); // fraction of ion density change in collapse time
  double ni(0.0); // ion number density
  double niz(0.0); // impurity ion number density
  double Zi(18.0); // impurity Z
  double Zeff(1.0); // effective Z
  double gamma_c(0.3); // thermal collapse rate (inverse timescale)
  double gamma_g(0.3); // solution dependent source rise rate 
  double gamma_src(0.3); // source rise rate
  double gmax(0.0); // solution dependent source factor
  double srcmax(0.0); // source max
  double srcsig(0.01); // source sigma in exponential (width)
  double t0(0); // time to turn on the temperature and density collapse
  double tf(0); // time to turn off the source
  double t(0.0); // time
  double ne(0.0); // electron number density
  double nem(0.0); // electron number density from previous time
  double ne0(0.0); // initial electron number density
  double jpar(0.0); // parallel current
  double jpartm(0.0); // parallel current from the previous timestep
  int jcon = 0; // conserve current after t0
  int setcurrent = 1; // current record switch
  int ionfunc = 0; // the time advance method for the ions

  // ---------------------------------------------------------------
  // Define the computational parameters 
  // ---------------------------------------------------------------

  int nt(200); // number of iterations
  int pltout(0); // realtime plotting option
  
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
      if (type == "nt") {
        in>>nt;
      } else if (type == "dtau") {
        in>>val;
        dtau=val;
      } else if (type == "pltout") {
        in>>pltout;
      } else if (type == "srcmax") {
        in>>srcmax;
      } else if (type == "srcsig") {
        in>>srcsig;
      } else if (type == "gamma_src") {
        in>>gamma_src;
      } else if (type == "gmax") {
        in>>gmax;
      } else if (type == "gamma_g") {
        in>>gamma_g;
      } else if (type == "gamma_c") {
        in>>gamma_c;
      } else if (type == "Tnorm") {
        in>>Tnorm;
      } else if (type == "Ti0") {
        in>>Ti0;
      } else if (type == "Tif") {
        in>>Tif;
      } else if (type == "ni0") {
        in>>ni0;
      } else if (type == "nif") {
        in>>nif;
      } else if (type == "nu") {
        in>>val;
        nu=val;
      } else if (type == "nui") {
        in>>val;
        nui=val;
      } else if (type == "Efield") {
        in>>val;
        Efield=val;
      } else if (type == "t0") {
        in>>t0;
      } else if (type == "tf") {
        in>>tf;
      } else if (type == "jcon") {
        in>>jcon;
      } else if (type == "ionfunc") {
        in>>ionfunc;
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
  Forms::Form_weightedIntegral current_form(mesh,fprev,Vpar,Jac); // form for computing jpar directly
  Forms::Form_weightedIntegral density_form(mesh,fprev,one,Jac); // form for computing ne directly

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
  a_kinetic.phiz = phiz;
  a_kinetic.psiz = psiz;
  a_kinetic.dtau = dtau;
  a_kinetic.E = Efield;
  a_kinetic.nu = nu;
  a_kinetic.nui = nui;
  a_kinetic.nuz = nuz;
  a_kinetic.mu = mu;
  a_kinetic.muz = muz;
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

  density_form.f=fprev;
  density_form.k=one;
  density_form.J=Jac;

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
  // the redundant variables should be collapsed into one
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
  phii_state.t0=t0;
  psii_state.t0=t0;
  psii_state.ionfunc=ionfunc;
  phii_state.ionfunc=ionfunc;
  phii_state.compute_coeffs();
  psii_state.compute_coeffs();
  phii.interpolate(phii_state); 
  psii.interpolate(psii_state); 
  
  // ---------------------------------------------------------------------
  // Set impurity ion potential parameters
  // the redundant variables should be collapsed into one
  // ---------------------------------------------------------------------
  phiz_state.phi0=niz;
  psiz_state.psi0=niz;
  phiz_state.phif=nif;
  psiz_state.psif=nif;
  phiz_state.Ti0=Ti0; 
  psiz_state.Ti0=Ti0; 
  phiz_state.Tif=Tif;
  psiz_state.Tif=Tif;
  phiz_state.Tnorm=Tnorm;
  psiz_state.Tnorm=Tnorm;
  phiz_state.mu=muz;
  psiz_state.mu=muz;
  phiz_state.gamma_phi=gamma_c;
  psiz_state.gamma_psi=gamma_c;
  phiz_state.t0=t0;
  psiz_state.t0=t0;
  psiz_state.ionfunc=ionfunc;
  phiz_state.ionfunc=ionfunc;
  phiz_state.compute_coeffs();
  psiz_state.compute_coeffs();
  phiz.interpolate(phiz_state); 
  psiz.interpolate(psiz_state); 

  // ---------------------------------------------------------------------
  // Set up solution independent source
  // ---------------------------------------------------------------------
  s.gamma_src=gamma_src;
  s.srcmax=srcmax;
  s.srcsig=srcsig;
  s.t0=t0;
  s.tf=tf;
  s.compute_coeffs();

  // ---------------------------------------------------------------------
  // Set up solution dependent source
  // ---------------------------------------------------------------------
  gs.gamma_g=gamma_g;
  gs.gmax=gmax;
  gs.t0=t0;
  gs.compute_coeffs();
  g.interpolate(gs);

  // ---------------------------------------------------------------------
  // Project the initial state to the finite element space and take a copy
  // ---------------------------------------------------------------------
  f.interpolate(initial_state);
  fprev=f;

  // ---------------------------------------------------------------------
  // Calculate the initial density and save it
  // ---------------------------------------------------------------------
  ne=assemble(density_form);
  ne0=ne;
  ni0=ne; //This goes against the input ni0 earlier, sort this out
  ni=ni0;

  // ---------------------------------------------------------------------
  // Set up a file for storing the solution (for creating a movie with paraview)
  // The rename calls ensure the variables are named correctly in the files.
  // ---------------------------------------------------------------------

  File file_f("f.pvd");
  File file_phi("phi.pvd");
  File file_phii("phii.pvd");
  File file_phiz("phiz.pvd");
  File file_psi("psi.pvd");
  File file_psii("psii.pvd");
  File file_psiz("psiz.pvd");
  std::ofstream file_dis("discharge.dat");
  f.rename("f","f");
  phi.rename("phi","phi"); 
  psi.rename("psi","psi");
  phii.rename("phii","phii");
  psii.rename("psii","psii");
  phiz.rename("phiz","phiz");
  psiz.rename("psiz","psiz");

  // ---------------------------------------------------------------------
  // Loop over time, this should be extremely simple now because all 
  // functions are already assigned to the form.
  // ---------------------------------------------------------------------

  swatch.stop();
  swatch = Timer("Time loop and problem solving");
  
  for (int it = 1; it <= nt; it++) {

    std::cout<<"time step: "<<it<<"/"<<nt<<" time: "<<t<< std::endl;

    // ---------------------------------------------------------------------
    // On the fly plotting for simple studies, use paraview to generate 
    // movies.
    // ---------------------------------------------------------------------
    
    switch (pltout) {
      case 1: plot(f,std::string("distribution function f"));
           // break;
      case 2: plot(g,std::string("sink function operator"));
           // break;
      case 3: plot(phii,std::string("phii_state"));
           // break;
      case 4: plot(psii,std::string("psii_state"));
           // break;
      case 5: plot(phiz,std::string("phiz_state"));
           // break;
      case 6: plot(psiz,std::string("psiz_state"));
           // break;
    }
    //if (pltout!=0) interactive();

    // ---------------------------------------------------------------------
    // Calculate the electron number density
    // ---------------------------------------------------------------------
    nem=ne;
    ne=assemble(density_form);

    // ---------------------------------------------------------------------
    // Set ion potentials
    // ---------------------------------------------------------------------
    if (ionfunc==0) {
      ni=phii_state.phii0;
      niz=phiz_state.phiz0;
    } else {
      // this method calculates the ion impurity density niz
      // from the ne source, then sets the density coefficient in the ion
      // distributions.  It works by a rule, for increasing density
      // niz is increased, for decreasing, niz and ni0 are equally decreased.
      if (ne-nem>=0.0) {
        niz=niz+(ne-nem)/Zi;
      } else {
        ni=ni*ne/nem;
        niz=niz*ne/nem;
      }
      phii_state.phi0=ni;
      psii_state.psi0=ni;
      phiz_state.phi0=niz;
      psiz_state.psi0=niz;
    }
    phii_state.t=t;
    psii_state.t=t;
    phiz_state.t=t;
    psiz_state.t=t;
    phii_state.compute_coeffs();
    psii_state.compute_coeffs();
    phiz_state.compute_coeffs();
    psiz_state.compute_coeffs();
    phii.interpolate(phii_state); 
    psii.interpolate(psii_state); 
    phiz.interpolate(phiz_state); 
    psiz.interpolate(psiz_state); 
 
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
    if (t>=t0 && setcurrent==1) {
      jpartm=jpar;
      setcurrent = 0;
    }
    jpar=assemble(current_form);

    // ---------------------------------------------------------------------
    // After t0 set the electric field to conserve the parrallel current
    // ---------------------------------------------------------------------
    if (jcon==1 && t>=t0 && jpar!=0.0) {
      Efield = Efield*jpartm/jpar;  // Is this the best way to constrain?
                                    // think this through
    }
    a_kinetic.E = Efield;

    // ---------------------------------------------------------------------
    // Write time stamps of f, phi, psi etc. to separate files. 
    // ---------------------------------------------------------------------
    if (it % 10==0) {
      file_f << f, t;
      file_phi << phi, t;
      file_psi << psi, t;    
      file_phii << phii, t;
      file_psii << psii, t;    
      file_phiz << phiz, t;
      file_psiz << psiz, t;    
    }

    file_dis << it <<" "<< t <<" "<< dtau <<" "<< jpar <<" "<< Efield 
             <<" "<< phii_state.Ti <<" "<< ni <<" "<< niz <<" "<< ne << "\n";
	
    // ---------------------------------------------------------------------
    // Advance the time
    // ---------------------------------------------------------------------
    t= t+dtau;
  }
  swatch.stop();
  list_timings();
  
  file_dis.close();

  if (pltout == 1) interactive();
  return 0;
}


