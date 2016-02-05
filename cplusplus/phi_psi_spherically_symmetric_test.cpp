#include <dolfin.h>

using namespace dolfin;

// Sub domain for applying Dirichlet boundary condition
class VeloBoundary : public SubDomain
{
private:
  static const double tol=1.e-15;
  double rightbound;  // this is not generic, should find a better way
  
public:
  // function used in fenics to flag boundary elements
  bool inside(const Array<double>& x, bool on_boundary) const{
    return on_boundary and abs(x[0]-rightbound)<tol;
  }
  
  // constructor given the upper limit (again, not generic)
  explicit VeloBoundary(const double& ul):rightbound(ul){}
};

// Source term (right-hand side)
class Source : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  { values[0] = exp(-x[0]*x[0]); }
};

// Jacobian term
class Jacobian : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  { values[0] = x[0]*x[0]; }
};

// The logarithmically space mesh (check it because looks funny)
class LogarithmicIntervalMesh : public IntervalMesh
{
  explicit LogarithmicIntervalMesh(const size_t& dim,const double& leftb, const double& rightb) : IntervalMesh(dim,leftb,rightb)
  {
    double pivot(log10(rightb)/double(dim));
    for(int i=0;i<dim-1;i++){
      coordinates[i] = pow(10,-2. + double(i)*pivot); // check this
    }
  }
};

// main programm
int main()
{
  // parameters for mesh creation (could appear as argv, argc)
  size_t nx=30;
  double xmin=0;
  double xmax=10;
  
  // Create mesh
  LogarithmicIntervalMesh mesh(nx-1,xmin,xmax);
  
  // Define function space for the mesh using Continuous Galerkin
  // (Lagrange) functions of order p on each element
  SphericallySymmetric::FunctionSpace V(mesh);
  
  // boundary object
  VeloBoundary boundary(xmax);
  // boundary condition
  //constant is evaluated as TruncGreenIntegral(xmin,xmax)
  Constant u0(0.0);

  DirichletBC bc(V, u0, boundary);

  Source f;
  Jacobian x2;

  // Define variational problem
  SphericallySymmetric::BilinearForm a(V, V);
  SphericallySymmetric::LinearForm L(V);
  L.f = f;
  L.J = x2;
  a.J = x2;

  // Compute solution
  Function u(V);
  solve(a == L, u, bc);

  // Save solution in VTK format
  File file("biharmonic.pvd");
  file << u;

  // Plot solution
  plot(u);
  interactive();

  return 0;
}
