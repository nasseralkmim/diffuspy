// Gmsh project created on Mon Apr 27 20:45:34 2015
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Physical Line(5) = {4};
Physical Line(6) = {1};
Physical Line(7) = {2};
Physical Line(8) = {3};
Line Loop(9) = {4, 1, 2, 3};
Plane Surface(10) = {9};
Physical Surface(11) = {10};
