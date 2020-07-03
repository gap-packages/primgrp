# This tests searches for primitive groups which are non-isomorphic,
# and have different structure descriptions, but have the same name
#
# In this version of the test the search is limited up to degree 64.
#
gap> START_TEST("duplicates.tst");

#
# Disable warnings which depend on Conway Polynomial databases
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
gap> ReadPackage("primgrp", "lib/testutils.g");
true
gap> PrimGrpNamesDuplicates(1,64);
(A(5) x A(5)):2 : [ [ 25, 23 ], [ 36, 17 ] ]
(S(5) x S(5)):2 : [ [ 25, 26 ], [ 36, 20 ] ]
2^6:(3 wreath Alt(3)):2 : [ [ 64, 22 ], [ 64, 23 ] ]
  [ 64, 22 ] non-isomorphic to [ 64, 23 ] : 
    ((A4 x A4 x A4) : C2) : C3
    ((C2 x C2 x C2 x C2 x C2 x C2) : ((C3 x C3 x C3) : C3)) : C2
2^6:3^2:Sym(3) : [ [ 64, 9 ], [ 64, 10 ] ]
  [ 64, 9 ] non-isomorphic to [ 64, 10 ] : 
    ((C2 x C2 x C2 x C2 x C2 x C2) : ((C3 x C3) : C3)) : C2
    (((C2 x C2 x C2 x C2 x C2 x C2) : (C3 x C3)) : C2) : C3
2^6:3^3:Sym(4) : [ [ 64, 35 ], [ 64, 36 ] ]
7^2:Q(8):3 : [ [ 49, 17 ], [ 49, 18 ] ]
A(10) : [ [ 10, 8 ], [ 45, 6 ] ]
A(5) : [ [ 5, 4 ], [ 10, 1 ] ]
A(6) : [ [ 6, 3 ], [ 15, 2 ] ]
A(7) : [ [ 7, 6 ], [ 15, 1 ], [ 21, 2 ], [ 35, 3 ] ]
A(8) : [ [ 8, 6 ], [ 28, 7 ], [ 35, 1 ] ]
A(9) : [ [ 9, 10 ], [ 36, 11 ] ]
Alt(5) wreath Sym(2) : [ [ 60, 3 ], [ 60, 4 ] ]
M(10) : [ [ 10, 6 ], [ 36, 3 ], [ 45, 2 ] ]
M(11) : [ [ 11, 6 ], [ 12, 1 ], [ 55, 4 ] ]
PGL(2, 11) : [ [ 12, 4 ], [ 55, 2 ], [ 55, 3 ] ]
PGL(2, 7) : [ [ 8, 5 ], [ 21, 1 ], [ 28, 1 ] ]
PGL(2, 9) : [ [ 10, 4 ], [ 36, 4 ], [ 45, 1 ] ]
PGammaL(2, 8) : [ [ 9, 9 ], [ 28, 3 ], [ 36, 2 ] ]
PGammaL(2, 9) : [ [ 10, 7 ], [ 36, 5 ], [ 45, 3 ] ]
PGammaU(3, 3) : [ [ 28, 5 ], [ 36, 7 ] ]
PSL(2, 11) : [ [ 12, 3 ], [ 55, 1 ] ]
PSL(2, 8) : [ [ 9, 8 ], [ 28, 2 ], [ 36, 1 ] ]
PSU(3, 3) : [ [ 28, 4 ], [ 36, 6 ], [ 63, 1 ], [ 63, 3 ] ]
PSU(3, 3).2 : [ [ 63, 2 ], [ 63, 4 ] ]
PSp(4, 3) : [ [ 27, 12 ], [ 36, 8 ], [ 40, 1 ], [ 40, 3 ], [ 45, 4 ] ]
PSp(4, 3):2 : [ [ 27, 13 ], [ 36, 9 ], [ 40, 2 ], [ 40, 4 ], [ 45, 5 ] ]
PSp(6, 2) : [ [ 28, 6 ], [ 36, 10 ], [ 63, 5 ] ]
S(10) : [ [ 10, 9 ], [ 45, 7 ] ]
S(5) : [ [ 5, 5 ], [ 10, 2 ] ]
S(6) : [ [ 6, 4 ], [ 15, 3 ] ]
S(7) : [ [ 7, 7 ], [ 21, 3 ], [ 35, 4 ] ]
S(8) : [ [ 8, 7 ], [ 28, 8 ], [ 35, 2 ] ]
S(9) : [ [ 9, 11 ], [ 36, 12 ] ]

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "duplicates.tst", 1);
