"""

This script tests the implementation of
the spherically symmetric nonlinear
Fokker-Planck equation with a spherically
symmetric Gaussian source term. The script
can be used to investigate how the source
rate and the source width affect the evolution
of the distribution function.
 
"""

import matplotlib.pyplot as plt
import scipy.special as sp
import numpy as np
from dolfin import *

# Variables to define the mesh,
# time step, iteration parameter
# and polynomial space degree
nx = 100   # number of grid points
xmin=0    # grid start point
xmax=10   # grid end point
dt = 0.1 # time step*collision frequency
nt = 100   # number of time steps
poly_order = 3  # polynomial degree of the basis functions         

# Construct the mesh and polynomial space
x = np.append(xmin,np.logspace(-3,np.log10(xmax),nx-1)).reshape(nx,1) # logspaced
x[nx-1]=xmax # make sure the last point is correct
mesh = IntervalMesh(nx-1,xmin,xmax) # generate the mesh
mesh.coordinates()[:]=x # set up the coordinates
V = FunctionSpace(mesh,"CG",poly_order) # continuous Galerkin elements of order poly_order

# Define the Dirichlet boundary domain at x=xmax
def Dirichlet_boundary_xmax(x,on_boundary):
    return on_boundary and abs(x[0]-xmax)<1E-15

# Define the boundary for projecting the initial state to FEM basis
def Dirichlet_boundary_initial(x,on_boundary):
    return on_boundary and (abs(x[0]-xmax)<1E-15 or abs(x[0]-xmin)<1E-15)


# Describe the initial state as a C++ class and
# create a dolfing Expression using the string
initial_state_str='''
namespace dolfin {
    class fzero : public Expression
    {

       public :
       fzero () : Expression () {}

       void eval (Array <double >& values , const Array <double >& x) const
       {
           values [0] = exp ( -x[0] * x[0] ) ;
       }
    };
}
'''

# Project the initial state to the function space
f0_expression=Expression(initial_state_str)
f0=TrialFunction(V)
source = Expression("1000*0.4*31*exp(-x[0]*x[0]/(2*0.0005))")
v=TestFunction(V)
a=inner(f0,v)*dx
L=f0_expression*v*dx
bc_initial=DirichletBC(V,f0_expression,Dirichlet_boundary_initial)
f0=Function(V)
solve(a == L, f0, bc_initial)


# take a copy of the initial state that we use in the iteration
fprev = f0

# an expression we need ever where
xsq = Expression("x[0]*x[0]")

#print assemble(xsq*xsq*f0*dx)/assemble(xsq*f0*dx)

# Generate the plot of the time evolution, plot the initial state
# with a thick line
xs=x.flatten()
plt.plot(xs,f0.compute_vertex_values(),linewidth=3)

# initialize arrays for storing temperature and density
Temp=np.zeros(nt)
dens=np.zeros(nt)

# loop over different time instances, use only one Picard iteration
for it in xrange(1,nt+1):
    
    # Compute the values of the potentials at the boundary
    phic=-assemble(xsq*fprev*dx)/xmax
    psic=-xmax/2*assemble(xsq*fprev*dx)-assemble(xsq*xsq*fprev*dx)/(6*xmax)

    # Define and solve for the potential functions using the boundary values
    # phi:
    bc_phi = DirichletBC(V,phic,Dirichlet_boundary_xmax)
    phi = TrialFunction(V)
    aphi = inner(-xsq*grad(phi),grad(v))*dx
    Lphi = fprev*v*xsq*dx
    phi = Function(V)
    solve(aphi == Lphi, phi, bc_phi)

    # psi:
    bc_psi = DirichletBC(V,psic,Dirichlet_boundary_xmax)
    psi = TrialFunction(V)
    apsi = inner(-xsq*grad(psi),grad(v))*dx
    Lpsi = phi*v*xsq*dx
    psi = Function(V)
    solve(apsi == Lpsi, psi, bc_psi)

    # Define and solve the weak problem for f
    fc = 0.
    bc_f=DirichletBC(V,fc,Dirichlet_boundary_xmax)
    f = TrialFunction(V)
    af = inner(xsq*v,f)*dx + dt*inner(xsq*grad(v),grad(phi)*f)*dx - dt*inner(xsq*grad(v),grad(grad(psi))*grad(f))*dx
    Lf = xsq*v*(fprev+dt*source)*dx
    f=Function(V)
    solve(af == Lf, f, bc_f)

    # plot and update the distribution function
    fprev = f
    if it == nt :
        plt.plot(xs,f.compute_vertex_values(),linewidth=3,color='red')
    else :
        plt.plot(xs,f.compute_vertex_values(),linewidth=1)

    dens[it-1] = assemble(xsq*f*dx);
    Temp[it-1] = assemble(xsq*xsq*f*dx)/dens[it-1];
    

# plot the evolution of the distribution
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-10,1e2])
plt.xlim([1e-1,1e1])
plt.xlabel('v')
plt.ylabel('f')   
plt.show()

# plot the density and temperature
plt.plot(dens,label='density')
plt.plot(Temp,label='temperature')
plt.xlabel('time step')
plt.legend()
plt.show()
