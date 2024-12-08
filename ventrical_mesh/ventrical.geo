h = 0.6; //mesh step

a_endo = 17.0;
b_endo = 7.0;
a_epi = 20.0;
b_epi = 10.0;
a_trunc = 5.0;
a_mid = 18.5;
b_mid = 8.5;

b_endo_trunc = b_endo*Sqrt(1.0 - a_trunc^2/a_endo^2);
b_epi_trunc  = b_epi*Sqrt(1.0 - a_trunc^2/a_epi^2);
b_mid_trunc  = b_mid*Sqrt(1.0 - a_trunc^2/a_mid^2);

Point(1)  = { b_endo,  0.0,    0.0, h};
Point(2)  = { 0.0   ,  0.0, a_endo, h};
Point(3)  = {-b_endo,  0.0,    0.0, h};
Point(4)  = { 0.0   ,  0.0,-a_endo, h};

Point(5)  = { b_epi,   0.0,    0.0, h};
Point(6)  = { 0.0   ,  0.0,  a_epi, h};
Point(7)  = {-b_epi,   0.0,    0.0, h};
Point(8)  = { 0.0   ,  0.0, -a_epi, h};

Point(10) = { 0.0,     0.0,    0.0, h};

Point(11) = {b_endo_trunc, 0.0, a_trunc, h};
Point(12) = {b_epi_trunc,  0.0, a_trunc, h};

Point(13)  = { -b_mid,   0.0,    0.0, h};
Point(14) =  {-b_mid_trunc,  0.0, a_trunc, h};
Point(15)  = { 0.0   ,  0.0,  -a_mid, h};

Ellipse(1) = {4, 10, 2, 1};
Ellipse(2) = {1, 10, 3, 11};
Line(3)    = {11, 12};
Ellipse(4) = {12, 10, 8, 5};
Ellipse(5) = {5, 10, 7, 8};
Line(6)    = {8, 15};
Ellipse(7) = {15, 10, 13, 14};
Line(8)    = {15, 4};

Line Loop(9)  = {1, 2, 3, 4, 5, 6, 8};
Plane Surface(10) = {9};

vol01[] = Extrude{{0, 0, 1}, {0, 0, 0}, 2*Pi/3}{ Surface{10}; };
vol02[] = Extrude{{0, 0, 1}, {0, 0, 0}, 2*Pi/3}{ Surface{vol01[0]}; };
vol03[] = Extrude{{0, 0, 1}, {0, 0, 0}, 2*Pi/3}{ Surface{vol02[0]}; };

Point {14} In Surface{vol02[4]};
Curve {7} In Volume{vol02[1]};

//command Physical create marks of surfaces, points, a etc that are avaliable in INMOST from tag GMSH_TAGS   
//1 - Dirichlet, 2 - endo, 3 - epi
Physical Surface ("endo", 2) = {vol01[2], vol01[3], vol02[2], vol02[3], vol03[2], vol03[3]};
Physical Surface ("base", 1) = {vol01[4], vol02[4], vol03[4]};
Physical Surface ("epi", 3) = {vol01[5], vol01[6], vol02[5], vol02[6], vol03[5], vol03[6]};
Physical Point (2) = {4};
Physical Point (3) = {8};
//1 - watch-curve
Physical Curve("midline", 1) = {7};


Physical Volume(99) = {vol01[1], vol02[1], vol03[1]};




