# primgrp, chapter 1
gap> START_TEST( "primgrp01.tst");

# primgrp/doc/../lib/primitiv.gd:138-149
gap> NrPrimitiveGroups(25);
28
gap> PrimitiveGroup(25,19);
5^2:((Q(8):3)'4)
gap> PrimitiveGroup(25,20);
ASL(2, 5)
gap> PrimitiveGroup(25,22);
AGL(2, 5)
gap> PrimitiveGroup(25,23);
(A(5) x A(5)):2

# primgrp/doc/prim.xml:34-42
gap> COHORTS_PRIMITIVE_GROUPS[49];
[ [ rec( parameter := 7, series := "Z", width := 2 ), 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
          20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33 ] ], 
  [ rec( parameter := [ 2, 7 ], series := "L", width := 2 ), [ 34 ] ], 
  [ rec( parameter := 7, series := "A", width := 2 ), [ 35, 36, 37, 38 ] ], 
  [ rec( parameter := 49, series := "A", width := 1 ), [ 39, 40 ] ] ]

# primgrp/doc/../lib/primitiv.gd:315-318
gap> PrimitiveIdentification(Group((1,2),(1,2,3)));
2

# primgrp/doc/../lib/primitiv.gd:249-254
gap> g:=PrimitiveGroup(25,2);
5^2:S(3)
gap> SimsNo(g);
3
gap> STOP_TEST("primgrp01.tst", 1 );
