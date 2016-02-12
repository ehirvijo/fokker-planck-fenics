lc1=2.0; // characteristic size at the half-circle
lc2=0.05; // characteristic size at the center
r = 10; // half circle radius

// Define four points from which the halfcircle is constructed
Point(1) = {0.0, 0.0, 0.0, lc2};
Point(2) = {0.0, -r, 0.0, lc1};
Point(3) = {r, 0.0, 0.0 ,lc1};
Point(4) = {0.0, r, 0.0, lc1};

// construct lines and arcs from the points
Line(1) = {1,2};
Circle(2) = {2,1,3};
Circle(3) = {3,1,4};
Line(4) = {4,1};

// construct the perimeter from the lines and the arcs
// to define the plane surface, i.e. the computational domain
Line Loop(5) = {1,2,3,4};

// construct the surface mesh to be inside the perimeter
Plane Surface(6) = {5};

// define the physical boundaries that dolfin should recognize
Physical Surface(0) = {6}; // the computational domain
Physical Line(0) = {4,1};  // the straight line of the half circle
Physical Line(1) = {2,3};  // the arc of the half circle

