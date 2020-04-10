//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {4, 0, 0, 1.0};
//+
Point(3) = {4, 4, 0, 1.0};
//+
Point(4) = {0, 4, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Physical Curve("bottom") = {1};
//+
Physical Curve("right") = {2};
//+
Physical Curve("top") = {3};
//+
Physical Curve("left") = {4};
//+
Surface(1) = {1};
//+
Physical Surface("fluid") = {1};
Transfinite Line "*" = 5 Using Bump 1.0;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
//+
Physical Surface("fluid") += {1};
