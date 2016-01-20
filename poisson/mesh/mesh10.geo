// Gmsh project created on Fri Feb 13 20:14:54 2015
Point(1) = {0, 0, -0, 1.0};
Point(2) = {10, 0, -0, 1.0};
Point(3) = {10, 10, -0, 1.0};
Point(4) = {0, 10, -0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Physical Line(7) = {4};
Physical Line(8) = {1};
Physical Line(9) = {2};
Physical Line(10) = {3};
Physical Surface(11) = {6};
