""" 

Solve the Rosenbluth potentials u1(x) and u2(x) 
for a Maxwellian plasma f(x)=exp(-x^2), i.e., solve

(1/x^2)*(x^2*u1')' = f(x) 
(1/x^2)*(x^2*u2')' = u1(x) 

which have the analytical solutions

u1(x) = -sqrt(pi)*erf(x)/(4*x)
u2(x) = -(2*x*exp(-x^2)-sqrt(pi)*(1+2*x^2)*erf(x))/(16*x)

Numerical solution requires truncation of the domain
so we use the truncated Green's function solutions
at the outer boundary x=c, resulting in conditions

u1'(0) = 0
u1(c)  = -(1/c)*Nintegrate[t^2*f(t),{t,0,c}]

u2'(0) = 0
u2(c)  = -(c/2)*Nintegrate[t^2*f(t),{t,0,c}]-(1/6c)*Nintegrate[t^4*f(t),{t,0,c}]

The outer boundary condition is computed with a
numerical quadrature to mimic a case where f would
be given as a polynomial.

"""

import matplotlib.pyplot as plt
import scipy.special as sp
import numpy as np
from dolfin import *
from scipy.integrate import quad

# number of mesh points, xmin and xmax
nx = 30   
xmin=0
xmax=10

# the actual mesh points, logarithmically spaced
x = np.append(xmin,np.logspace(-2,np.log10(xmax),nx-1)).reshape(nx,1)
x[nx-1]=xmax

# Polynomial order of trial/test functions
p = 3           

# create mesh 
mesh = IntervalMesh(nx-1,xmin,xmax)
mesh.coordinates()[:]=x

# Define function space for the mesh using Continuous Galerkin
# (Lagrange) functions of order p on each element
V = FunctionSpace(mesh,"CG",p)

# Compute the truncated Green's function solutions at the boundary
u1_c = -quad(lambda t: t**2*np.exp(-t**2),xmin,xmax)[0]/xmax
u2_c = -xmax/2*quad(lambda t: t**2*np.exp(-t**2),xmin,xmax)[0]-quad(lambda t: t**4*np.exp(-t**2),xmin,xmax)[0]/(6*xmax)

# Set up the "Dirichlet boundary domain"
def Dirichlet_boundary(x,on_boundary):
    return on_boundary and abs(x[0]-xmax)<1E-15

# SET UP AND SOLVE THE PROBLEM 1
bc_u1 = DirichletBC(V,u1_c,Dirichlet_boundary)
u1 = TrialFunction(V)
v = TestFunction(V)
x2 = Expression("x[0]*x[0]")
f = Expression("exp(-x[0]*x[0])")
a1 = inner(-x2*grad(u1),grad(v))*dx
L1 = f*v*x2*dx

u1 = Function(V)
solve(a1 == L1, u1, bc_u1)

## SET UP AND SOLVE THE PROBLEM 2
bc_u2 = DirichletBC(V,u2_c,Dirichlet_boundary)
u2 = TrialFunction(V)
v = TestFunction(V)
x2 = Expression("x[0]*x[0]")
a2 = inner(-x2*grad(u2),grad(v))*dx
L2 = u1*v*x2*dx

u2 = Function(V)
solve(a2 == L2, u2, bc_u2)

# evaluate the solutions at mesh vertices and compare
# with the analytical expressions
xs=x.flatten()
plt.plot(xs,-np.sqrt(np.pi)*sp.erf(xs)/(4*xs),linewidth=2,label='Analytical phi')
plt.plot(xs,-(2*np.exp(-xs**2)+np.sqrt(np.pi)*(1.+2*xs**2)*sp.erf(xs)/xs)/16,linewidth=2,label='Analytical psi')
plt.plot(xs,u1.compute_vertex_values(),linewidth=2,label='FEM phi')
plt.plot(xs,u2.compute_vertex_values(),linewidth=2,label='FEM psi')

#plt.plot(xs,grad(u2).compute_vertex_values(),linewidth=2,label='FEM psi')
plt.legend(loc=3)
plt.xscale('log')
plt.xlabel('v')
plt.ylabel('phi, psi')
plt.show()


