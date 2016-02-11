lc=0.3; // element "size"
r = 10; // circle radius

// Define four points from which the halfcircle is constructed
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {0.0, -r, 0.0, lc};
Point(3) = {r, 0.0, 0.0 ,lc};
Point(4) = {0.0, r, 0.0, lc};

// construct lines and arcs from the points
Line(1) = {1,2};
Ellipse(2) = {2,1,1,3};
Ellipse(3) = {3,1,1,4};
Line(4) = {4,1};

// construct the perimeter from the lines and the arcs
Line Loop(5) = {1,2,3,4};

// construct the surface mesh to be inside the perimeter
Plane Surface(6) = {5};