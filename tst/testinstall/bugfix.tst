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

# The name of this group was written "... 4 \circ 2^(1+4) ..." in gps4.g, and
# GAP reads \c as character 3, so the name came out as "4 <chr 3>irc". The
# subgroup really is the central product of C4 with the extraspecial group
# 2^(1+4), for which ATLAS notation is the circle.
gap> Name(PrimitiveGroup(625,657));
"5^4:(4 \\circ 2^(1+4)).Sp(4, 2)"
gap> ForAny([2..4095], function(d)
>      PrimGrpLoad(d);
>      return ForAny(PRIMGRP[d], e -> ForAny(e[7], c -> INT_CHAR(c) < 32));
>    end);
false

# The socle series for unitary groups is "2A"; 51 entries in degrees >= 2500
# spelled it "^2A", which GAP never produces.
gap> SocleTypePrimitiveGroup(PrimitiveGroup(2752,1));
rec( parameter := [ 3, 7 ], series := "2A", width := 1 )
gap> SocleTypePrimitiveGroup(PrimitiveGroup(2500,18));
rec( parameter := [ 2, 5 ], series := "2A", width := 2 )
gap> ForAny([2..4095], function(d)
>      PrimGrpLoad(d);
>      return ForAny(PRIMGRP[d], e -> e[8][1][1] = '^');
>    end);
false

# Twenty-nine socle types recorded a rank one too high, in degrees 3159 to
# 4095: PrimitiveGroup(4095,1) claimed C(7,2) where the socle is B(6,2).  The
# series letter moves with some of them, which is harmless where the two names
# denote the same group, but the rank never was.
gap> List([[3640,1],[3969,9],[4095,1]], p -> PRIMGrp(p[1],p[2])[8]);
[ [ "C", [ 3, 3 ], 1 ], [ "B", [ 3, 2 ], 2 ], [ "B", [ 6, 2 ], 1 ] ]

# one from each affected degree: what is stored must be what GAP computes
gap> soc := p -> SocleTypePrimitiveGroup(PrimitiveGroup(p[1],p[2]));;
gap> flat := s -> [s.series, s.parameter, s.width];;
gap> ForAll([[3159,3],[3240,9],[3280,1],[3321,14],
>            [3640,2],[3906,3],[3969,9],[4095,1]],
>           p -> flat(soc(p)) = PRIMGrp(p[1],p[2])[8]);
true

#
gap> STOP_TEST("bugfix.tst", 1);
