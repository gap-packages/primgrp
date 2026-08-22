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
gap> ForAny([2..4095], d -> ForAny([1..NrPrimitiveGroups(d)],
>            i -> ForAny(PRIMGrp(d,i)[7], c -> INT_CHAR(c) < 32)));
false

# The socle series for unitary groups is "2A"; 51 entries in degrees >= 2500
# spelled it "^2A", which GAP never produces.
gap> SocleTypePrimitiveGroup(PrimitiveGroup(2752,1));
rec( parameter := [ 3, 7 ], series := "2A", width := 1 )
gap> SocleTypePrimitiveGroup(PrimitiveGroup(2500,18));
rec( parameter := [ 2, 5 ], series := "2A", width := 2 )
gap> ForAny([2..4095], d -> ForAny([1..NrPrimitiveGroups(d)],
>            i -> PRIMGrp(d,i)[8][1][1] = '^'));
false

# Three entries claimed the socle was B(4,3) = O(9,3), of order
# 65784756654489600, when the group itself has order 4585351680 = |O(7,3)|.
# Found by the consistency check in tst/testutils.g.
gap> List([[3159,3],[3640,2],[3640,4]],
>         p -> SocleTypePrimitiveGroup(PrimitiveGroup(p[1],p[2])).parameter);
[ [ 3, 3 ], [ 3, 3 ], [ 3, 3 ] ]

# Sixteen more socle types with the rank one too high: C(3,9) where the socle
# has order 1721606400 = |O(5,9)| = |B(2,9)|, and B(5,3) where it has order
# 65784756654489600 = |O(9,3)| = |B(4,3)|.
gap> Set([[3240,9],[3240,13],[3321,9],[3321,13]],
>         p -> SocleTypePrimitiveGroup(PrimitiveGroup(p[1],p[2])));
[ rec( parameter := [ 2, 9 ], series := "B", width := 1 ) ]
gap> Set([[3240,14],[3240,15],[3280,2],[3280,4],[3321,14],[3321,15]],
>         p -> SocleTypePrimitiveGroup(PrimitiveGroup(p[1],p[2])));
[ rec( parameter := [ 4, 3 ], series := "B", width := 1 ) ]

# B(4,3) and C(4,3) have the same order, 65784756654489600, but are not
# isomorphic, so these two had to be told apart by series and not by order.
gap> List([[3280,1],[3280,2]],
>         p -> SocleTypePrimitiveGroup(PrimitiveGroup(p[1],p[2])));
[ rec( parameter := [ 4, 3 ], series := "C", width := 1 ), 
  rec( parameter := [ 4, 3 ], series := "B", width := 1 ) ]

# Degree 3640 carries both series at rank 3, of equal order 4585351680, and
# both were stored one rank too high.
gap> List([[3640,1],[3640,2],[3640,3],[3640,4]],
>         p -> SocleTypePrimitiveGroup(PrimitiveGroup(p[1],p[2])));
[ rec( parameter := [ 3, 3 ], series := "C", width := 1 ), 
  rec( parameter := [ 3, 3 ], series := "B", width := 1 ), 
  rec( parameter := [ 3, 3 ], series := "C", width := 1 ), 
  rec( parameter := [ 3, 3 ], series := "B", width := 1 ) ]

# The last six, including PrimitiveGroup(4095,1), which is the example that
# prompted issue #40: the socle was recorded as C(7,2) rather than B(6,2).
gap> List([[3906,3],[3906,4],[3969,9],[4095,1]],
>         p -> SocleTypePrimitiveGroup(PrimitiveGroup(p[1],p[2])));
[ rec( parameter := [ 3, 5 ], series := "C", width := 1 ), 
  rec( parameter := [ 3, 5 ], series := "B", width := 1 ), 
  rec( parameter := [ 3, 2 ], series := "B", width := 2 ), 
  rec( parameter := [ 6, 2 ], series := "B", width := 1 ) ]

# Type 4c entries stored as ["pa",m,b,topgens]: the top group's generators are
# permutations, and a generator that is a product of two or more cycles was
# once written out as several separate generators, giving a larger group.
gap> List([[625,6],[625,41],[1296,26]],
>         p -> Collected(List(OrbitsDomain(Stabilizer(PrimitiveGroup(p[1],p[2]),1),
>                                          [1..p[1]]), Length)));
[ [ [ 1, 1 ], [ 16, 1 ], [ 32, 3 ], [ 256, 2 ] ], 
  [ [ 1, 1 ], [ 16, 1 ], [ 32, 3 ], [ 256, 2 ] ], 
  [ [ 1, 1 ], [ 20, 1 ], [ 50, 3 ], [ 500, 1 ], [ 625, 1 ] ] ]

# Type 2 entries stored as ["cl",<series>,<d>,<q>,<ext>]: the socle acting on
# its natural geometry, extended by <ext> outer automorphisms.  Degree 126 is
# the case that shows why the point set is built by hand -- GAP's own
# PSU(3,5) acts on all 651 projective points, not on the 126 isotropic ones.
gap> g := PrimitiveGroup(126,3);;
gap> Size(g) = Size(PSU(3,5)) and NrMovedPoints(g) = 126 and IsPrimitive(g);
true
gap> List([[126,6],[364,9],[820,7],[1066,1]],
>         p -> Size(PrimitiveGroup(p[1],p[2])));
[ 756000, 42064805779476480, 6886425600, 10151968619520 ]

# Type 2 entries stored as ["cl"/"clw",...,<variety>].  B(2,q) acts on the
# (q^2+1)(q+1) singular points and on just as many totally singular lines,
# and for odd q those actions are inequivalent -- the quadrangle Q(4,q) is
# self-dual only in even characteristic.  So degree 40 carries PSp(4,3)
# twice, with equal order, equal socle and equal suborbits; what differs is
# the point stabiliser, the two maximal parabolics.
gap> List([[40,1],[40,3]], p -> Collected(List(OrbitsDomain(
>         Stabilizer(PrimitiveGroup(p[1],p[2]),1), [2..40]), Length)));
[ [ [ 12, 1 ], [ 27, 1 ] ], [ [ 12, 1 ], [ 27, 1 ] ] ]
gap> List([[40,1],[40,3]],
>         p -> AbelianInvariants(Stabilizer(PrimitiveGroup(p[1],p[2]), 1)));
[ [ 3 ], [ 2 ] ]

# Affine entries of prime degree are stored as PGPrime(d), one per divisor d
# of p-1.  Degrees below 4096 keep the names they always had; degrees 4096 to
# 8191 had no name at all before, and now use the same convention.
gap> List([[59,1],[59,2],[59,3],[59,4]], p -> PRIMGrp(p[1],p[2])[7]);
[ "C(59)", "D(2*59)", "59:29", "AGL(1, 59)" ]
gap> List([[4099,1],[4099,2],[4099,7]], p -> PRIMGrp(p[1],p[2])[7]);
[ "C(4099)", "D(2*4099)", "4099:2049" ]
gap> g := PrimitiveGroup(4099,7);;
gap> Size(g) = 4099*2049 and NrMovedPoints(g) = 4099 and IsPrimitive(g);
true

#
gap> STOP_TEST("bugfix.tst", 1);
