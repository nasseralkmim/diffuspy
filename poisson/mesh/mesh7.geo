// Gmsh project created on Fri Feb 13 09:58:10 2015
Point(1) = {0, 1, 0, 1.0};
Point(2) = {-0, -0, 0, 1.0};
Point(3) = {1.7, -0, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Physical Line(4) = {1};
Physical Line(5) = {2};
Physical Line(6) = {3};
Line Loop(7) = {3, 1, 2};
Plane Surface(8) = {7};
Physical Surface(9) = {8};
