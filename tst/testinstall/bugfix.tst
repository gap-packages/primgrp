gap> START_TEST("bugfix.tst");

# Use consistent names for the Tits group and for the Ree group containing it.
# See <https://github.com/gap-packages/primgrp/issues/17>
gap> G := PrimitiveGroup(1600, 20);
2F(4, 2)'
gap> IsomorphismTypeInfoFiniteSimpleGroup(G);
rec( name := "2F(4,2)' = Ree(2)' = Tits", parameter := 2, series := "2F", 
  shortname := "2F4(2)'" )

#
gap> G := PrimitiveGroup(1755, 1);
2F(4, 2)'
gap> IsomorphismTypeInfoFiniteSimpleGroup(G);
rec( name := "2F(4,2)' = Ree(2)' = Tits", parameter := 2, series := "2F", 
  shortname := "2F4(2)'" )

#
gap> G := PrimitiveGroup(1755, 2);
2F(4, 2)

#
gap> G := PrimitiveGroup(2304, 9);
2F(4, 2)'
gap> IsomorphismTypeInfoFiniteSimpleGroup(G);
rec( name := "2F(4,2)' = Ree(2)' = Tits", parameter := 2, series := "2F", 
  shortname := "2F4(2)'" )

#
gap> G := PrimitiveGroup(2304, 10);
2F(4, 2)

#
gap> G := PrimitiveGroup(625,672);
5^4:4.(S5 wr 2)

# Fix of the bug reported by Ali Reza Rahimipour on 12 June 2020
gap> G := PrimitiveGroup(1600,17);
PSU(3, 4)
gap> G := PrimitiveGroup(1600,18);
PSU(3, 4):2
gap> G := PrimitiveGroup(1600,19);
PSU(3, 4):4

# not really a bug, but: ensure consisting naming for O+(8,2) and friends;
# PrimitiveGroup(960,7) used to have the name "OPlus(8, 2)"
gap> degs := [ 119, 120, 135, 136, 765, 960, 1071, 1120, 1575, 1632 ];;
gap> AllPrimitiveGroups(NrMovedPoints, degs, Size, Size(SO(+1,8,2)));
[ PSO+(8, 2), PSO+(8, 2), PSO+(8, 2), POmega+(8, 2):2, POmega+(8, 2):2 ]
gap> AllPrimitiveGroups(NrMovedPoints, degs, Size, Size(SO(+1,8,2))/2);
[ O+(8, 2), O+(8, 2), O+(8, 2), POmega+(8, 2), POmega+(8, 2) ]
gap> AllPrimitiveGroups(NrMovedPoints, degs, Size, Size(SO(-1,8,2)));
[ PSO-(8, 2), PSO-(8, 2), PSO-(8, 2), PSO-(8,2), POmega-(8, 2):2 ]
gap> AllPrimitiveGroups(NrMovedPoints, degs, Size, Size(SO(-1,8,2))/2);
[ O-(8, 2), O-(8, 2), O-(8, 2), POmega-(8, 2), POmega-(8, 2) ]

# Fix an error in the name reported by Andries Brouwer
gap> G := PrimitiveGroup(81,103);
3^4:(2^3:Sym(4)):Sym(3)

#
gap> G := PrimitiveGroup(100,14);
Alt(6)^2.4

#
gap> G := PrimitiveGroup(121,38);
11^2:(5 x Q_24)

#
gap> G := PrimitiveGroup(289,35);
17^2:SL(2,3)

#
gap> STOP_TEST("bugfix.tst", 1);
