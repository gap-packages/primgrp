#############################################################################
##
##  performance.tst
##
##  Constructing a group must not fall back on searching a list of points.
##
##  The affine entries of degree 4096 and up enumerate F_p^d by the integer
##  the base-p digits spell out.  Looking each image up in that list is a
##  linear scan: one generator of PrimitiveGroup(6561,1) took 4985 ms that
##  way against 17 ms when the image number is computed from the image
##  vector.  Nothing else in the suite would notice, because this point
##  enumeration only occurs above degree 4095.
##
gap> START_TEST("performance.tst");

# The arithmetic action agrees with the one found by searching.  Done on a
# small space, so that the slow side of the comparison stays cheap.
gap> v := PRIMGRP_AffineVectors(5, 3, "int");;
gap> Length(v);
125
gap> m := ImmutableMatrix(GF(5), [[1,1,0],[0,1,1],[1,0,1]] * Z(5)^0);;
gap> PRIMGRP_AffineAction(5, 3, "int", v, m) = Permutation(m, v, OnRight);
true
gap> t := First(v, i -> not IsZero(i));;
gap> PRIMGRP_AffineTranslation(5, 3, "int", v, t)
>      = Permutation(t, v, function(a,b) return a+b; end);
true

# A generous ceiling on the degrees that actually use it.  Not a benchmark,
# just a tripwire: these take tens of milliseconds today.
gap> start := Runtime();;
gap> G := PrimitiveGroup(6561,1);;
gap> H := PrimitiveGroup(4096,1);;
gap> K := PrimitiveGroup(8191,1);;
gap> Runtime() - start < 30000;
true
gap> [NrMovedPoints(G), NrMovedPoints(H), NrMovedPoints(K)];
[ 6561, 4096, 8191 ]

#
gap> STOP_TEST("performance.tst", 1);
