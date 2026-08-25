gap> START_TEST("product4c.tst");

#
gap> ReadPackage("primgrp", "tst/product4c.g");
true

# A group of type 4c taken apart over its coordinates and put back together.
# The check is an identity -- every generator rebuilds as itself relabelled --
# so it holds whatever generators the random search happened to pick.
gap> List([[25,24],[36,13],[64,65]],
>         function(p)
>           local g, r;
>           g := PrimitiveGroup(p[1],p[2]);
>           r := PRIMGRP_DecomposeProductAction(g);
>           return [ r.m, r.k, PRIMGRP_CheckProductAction(g, r) ];
>         end);
[ [ 5, 2, true ], [ 6, 2, true ], [ 8, 2, true ] ]

# The number of coordinates may be given, for a group that is not an entry of
# the library and so has no socle type to read it off.
gap> g := Group(GeneratorsOfGroup(PrimitiveGroup(25,24)));;
gap> r := PRIMGRP_DecomposeProductAction(g, 2);;
gap> PRIMGRP_CheckProductAction(g, r);
true

#
gap> STOP_TEST("product4c.tst", 1);
