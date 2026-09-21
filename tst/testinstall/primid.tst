gap> START_TEST("primid.tst");
gap> ForAll([2..1000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([1001..2000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([2001..3000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([3001..4000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([4001..4095], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true

# Recognising Alt(n) or Sym(n) records what it found.
gap> x := PermList(Reversed([1..20]));;
gap> G := Group(List(GeneratorsOfGroup(SymmetricGroup(20)), y -> y^x));;
gap> PrimitiveIdentification(G);
4
gap> HasIsNaturalSymmetricGroup(G) and IsNaturalSymmetricGroup(G) and HasSize(G);
true
gap> G := Group(List(GeneratorsOfGroup(AlternatingGroup(20)), y -> y^x));;
gap> PrimitiveIdentification(G);
3
gap> HasIsNaturalAlternatingGroup(G) and IsNaturalAlternatingGroup(G) and HasSize(G);
true

# A library group knows its own number, so the tests above compute nothing.
# Conjugated generators give a group that has to be identified for real.
gap> roundtrip := function(n)
>   local x, i, g, h;
>   x := PermList(Reversed([1..n]));
>   for i in [1..NrPrimitiveGroups(n)] do
>     g := PrimitiveGroup(n, i);
>     h := Group(List(GeneratorsOfGroup(g), y -> y^x));
>     if PrimitiveIdentification(h) <> i then
>       return false;
>     fi;
>   od;
>   return true;
> end;;
gap> Filtered([2..100], n -> not roundtrip(n));
[  ]

# Affine groups that agree in every invariant tried: the point stabilisers
# are two embeddings of SmallGroup(208,3) into GL(4,5), not conjugate there.
gap> x := PermList(Reversed([1..625]));;
gap> List([528, 529], i -> PrimitiveIdentification(
>      Group(List(GeneratorsOfGroup(PrimitiveGroup(625, i)), y -> y^x))));
[ 528, 529 ]

# Product action.  The first seven agree in order, suborbits and transitivity
# and have socle Alt(6)^3 on 10^3 points; the next two have socle Alt(8)^4, and
# telling them apart by conjugacy classes took a quarter of an hour each.
gap> ident := function(n, i)
>   local x;
>   x := PermList(Reversed([1..n]));
>   return PrimitiveIdentification(
>     Group(List(GeneratorsOfGroup(PrimitiveGroup(n, i)), y -> y^x)));
> end;;
gap> List([88..94], i -> ident(1000, i));
[ 88, 89, 90, 91, 92, 93, 94 ]
gap> List([1448, 1454], i -> ident(4096, i));
[ 1448, 1454 ]

# PSL(3,5) on points and on lines: the socle factors of the group and of the
# library's copy need a conjugating permutation in Sym(31)
gap> PD := [ PRIMGrp(961, 1) ];;
gap> x := PermList(Reversed([1..961]));;
gap> G := Group(List(GeneratorsOfGroup(PrimitiveGroup(961, 1)), y -> y^x));;
gap> PRIMGRP_IdentifyProductAction(G, PD) = PD;
true

# no conjugating permutation: the groups are not isomorphic, or, for S4 on
# the cosets of C4 and of a non-normal 2^2, no isomorphism maps point
# stabilisers to point stabilisers
gap> PRIMGRP_ConjugatingPermutation(Group((1,2,3,4,5)), AlternatingGroup(5));
fail
gap> PRIMGRP_ConjugatingPermutation(TransitiveGroup(6, 7), TransitiveGroup(6, 8));
fail

# data/primid.g holds what PRIMGRP_IdClassData computes
gap> ForAll([9, 10, 16, 64, 81], deg -> PRIMGRP_IdData(deg) =
>      List(PRIMGRP_IdClasses(deg), c -> PRIMGRP_IdClassData(deg, c)));
true

# a candidate outside the stored class leaves the candidates unchanged
gap> PD := List([3, 1], i -> PRIMGrp(9, i));;
gap> PRIMGRP_FilterByCycleTypes(PrimitiveGroup(9, 3), PD) = PD;
true

# cycle counts that differ by the hash modulus collide
gap> PRIMGRP_IdClassDataFromTypes([1, 2], [1, 1], [[[1]], [[2147483648]]]);
Error, hash collision for groups [ 1, 2 ]
gap> STOP_TEST( "primid.tst", 1);
