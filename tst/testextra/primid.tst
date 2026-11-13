gap> START_TEST("primid.tst");

# data/primid.g holds what PRIMGRP_IdClassData computes; beyond degree 512
# that takes hours
gap> Filtered([2..512], deg -> PRIMGRP_IdData(deg) <>
>      List(PRIMGRP_IdClasses(deg), c -> PRIMGRP_IdClassData(deg, c)));
[  ]

# Size of D(2*8191) by Schreier-Sims takes minutes
gap> x := PermList(Reversed([1..8191]));;
gap> PrimitiveIdentification(
>      Group(List(GeneratorsOfGroup(PrimitiveGroup(8191, 2)), y -> y^x)));
2

# a random search for Jordan elements misses about once in a thousand
# conjugates of Sym(7526) in IsNaturalSymmetricGroup
gap> x := PermList(Reversed([1..7526]));;
gap> PrimitiveIdentification(
>      Group(List(GeneratorsOfGroup(SymmetricGroup(7526)), y -> y^x)));
2

# with socle Alt(73)^2: Size on the 5329 points takes minutes
gap> x := PermList(Reversed([1..5329]));;
gap> PrimitiveIdentification(
>      Group(List(GeneratorsOfGroup(PrimitiveGroup(5329, 280)), y -> y^x)));
280
gap> STOP_TEST("primid.tst", 1);
