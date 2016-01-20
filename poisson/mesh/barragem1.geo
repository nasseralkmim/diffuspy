// Gmsh project created on Tue Jul 14 13:03:05 2015
Point(1) = {0, 0, 0, 1.0};
Point(2) = {95.8, 0, 0, 1.0};
Point(3) = {18.513, 99.06, 0, 1.0};
Point(4) = {14.701, 106.68, 0, 1.0};
Point(5) = {14.701, 121.92, 0, 1.0};
Point(6) = {4.953, 121.92, 0, 1.0};
Point(7) = {4.953, 99.06, 0, 1.0};
Translate {0, .3, 0} {
  Point{5};
}
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
Line Loop(8) = {2, 3, 4, 5, 6, 7, 1};
Plane Surface(9) = {8};
Physical Line(10) = {1};
Physical Line(11) = {2};
Physical Line(12) = {3};
Physical Line(13) = {4};
Physical Line(14) = {5};
Physical Line(15) = {6};
Physical Line(16) = {7};
Physical Surface(17) = {9};
Delete {
  Line{1};
}
Delete {
  Surface{9};
}
Delete {
  Line{1};
}
Point(8) = {-100, 0, 0, 1.0};
Point(9) = {-100, -70, 0, 1.0};
Point(10) = {220, -70, 0, 1.0};
Point(11) = {220, 0, 0, 1.0};
Line(18) = {1, 8};
Line(19) = {8, 9};
Line(20) = {9, 10};
Line(21) = {10, 11};
Line(22) = {11, 2};
Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
Coherence;
Line Loop(23) = {7, 18, 19, 20, 21, 22, 2, 3, 4, 5, 6};
Plane Surface(24) = {23};
