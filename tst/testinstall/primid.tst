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

# data/primid.g holds what PRIMGRP_IdClassData computes
gap> ForAll([9, 10, 16, 64, 81], deg -> PRIMGRP_IdData(deg) =
>      List(PRIMGRP_IdClasses(deg), c -> PRIMGRP_IdClassData(deg, c)));
true
gap> STOP_TEST( "primid.tst", 1);
