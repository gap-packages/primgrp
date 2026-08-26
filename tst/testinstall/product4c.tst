gap> START_TEST("product4c.tst");

#
gap> ReadPackage("primgrp", "tst/product4c.g");
true

# Test PRIMGRP_DecomposeProductAction for a few sample inputs
gap> List([[25,24],[36,13],[64,65]],
>         function(p)
>           local g, r;
>           g := PrimitiveGroup(p[1],p[2]);
>           r := PRIMGRP_DecomposeProductAction(g);
>           return [ r.m, r.k, PRIMGRP_CheckProductAction(g, r) ];
>         end);
[ [ 5, 2, true ], [ 6, 2, true ], [ 8, 2, true ] ]

# The number of coordinates may be given, for a group that is not an entry of
# If SocleTypePrimitiveGroup is not set, then one has to supply the second
# argument to PRIMGRP_DecomposeProductAction; so test that case as well.
gap> g := Group(GeneratorsOfGroup(PrimitiveGroup(25,24)));;
gap> r := PRIMGRP_DecomposeProductAction(g, 2);;
gap> PRIMGRP_CheckProductAction(g, r);
true

#
gap> STOP_TEST("product4c.tst", 1);
