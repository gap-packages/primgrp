#############################################################################
##
##  performance.tst
##
##  A group is built once a session, not once per question asked about it.
##
##  PrimitiveGroup is called afresh for every attribute a caller wants, and
##  most of the library is now a construction rather than a list of
##  permutations.  PRIMGRP_GENCACHE keeps the generators of what was built, so
##  a second call assembles the group from those instead of running the
##  construction again: over degrees 2 to 500 a second pass cost 1154 ms and
##  costs 39 ms.
##
##  Nothing else in the suite would notice it stop working.  The rest of the
##  tests ask for a group once, or ask for the fields of an entry rather than
##  the group.
##
gap> START_TEST("performance.tst");

# A group from the cache is the group that was built: same generators, and the
# same answers to what the entry sets.
gap> G := PrimitiveGroup(2401,1);; H := PrimitiveGroup(2401,1);;
gap> [ GeneratorsOfGroup(G) = GeneratorsOfGroup(H), Size(G) = Size(H),
>      Transitivity(G) = Transitivity(H), ONanScottType(G) = ONanScottType(H) ];
[ true, true, true, true ]

# The natural alternating and symmetric groups are left out of the cache, and
# this is why: rebuilt from bare generators they would no longer know what they
# are, and finding out again costs more than building them did.
gap> G := PrimitiveGroup(20,3);; H := PrimitiveGroup(20,3);;
gap> [ Name(G), HasIsNaturalAlternatingGroup(G), HasIsNaturalAlternatingGroup(H) ];
[ "A(20)", true, true ]

# Twenty calls cost no more than the one that built the group.  A ratio rather
# than a wall clock, so it says the same thing on a slow machine.
gap> t := Runtime();; K := PrimitiveGroup(3125,1);; one := Runtime() - t;;
gap> t := Runtime();; for i in [1..20] do K := PrimitiveGroup(3125,1);; od;
gap> twenty := Runtime() - t;;
gap> twenty <= one;
true

#
gap> STOP_TEST("performance.tst", 1);
