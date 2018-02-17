gap> START_TEST("irrednumbers.tst");
gap> n := Filtered([2..255],IsPrimePowerInt);;
gap> n := List(n, Factors);;
gap> n := Filtered(n, t -> Length(t) > 1 );;
gap> n := List(n, t -> [ Length(t), t[1] ] );;
gap> List(n, t -> NumberIrreducibleSolvableGroups( t[1], t[2] ));
[ 2, 2, 7, 10, 19, 9, 2, 29, 40, 108, 42, 22, 2, 62, 16 ]
gap> Sum( List( n, t -> NumberIrreducibleSolvableGroups( t[1], t[2] )));
372
gap> ForAll(n, t -> NumberIrreducibleSolvableGroups( t[1], t[2] ) =
> Length( AllIrreducibleSolvableGroups( Dimension, t[1], Characteristic,t[2] )));
true
gap> ForAll(n, t -> IsSolvable(OneIrreducibleSolvableGroup(
>                                Dimension, t[1], Characteristic,t[2] )));
true
gap> STOP_TEST( "irrednumbers.tst", 1);
