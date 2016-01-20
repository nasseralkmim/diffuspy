// Gmsh project created on Thu Feb 05 00:23:47 2015
Point(1) = {-0.2, 0.9, -0, 1.0};
Point(2) = {-0.6, 0.3, -0, 1.0};
Point(3) = {0.7, 0.3, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line Loop(4) = {3, 1, 2};
Plane Surface(5) = {4};
Physical Line(6) = {1};
Physical Line(7) = {3};
Physical Line(8) = {2};
Physical Surface(9) = {5};
