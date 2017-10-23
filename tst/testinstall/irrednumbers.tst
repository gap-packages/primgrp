gap> START_TEST("irrednumbers.tst");
gap> n := Filtered([2..255],IsPrimePowerInt);;
gap> n := List(n, Factors);;
gap> n := Filtered(n, t -> Length(t) > 1 );;
gap> n := List(n, t -> [ Length(t), t[1] ] );;
gap> List(n, t -> NumberIrreducibleSolvableGroups( t[1], t[2] ));
[ 2, 2, 7, 10, 19, 9, 2, 29, 40, 108, 42, 22, 2, 62, 16 ]
gap> Sum( List( n, t -> NumberIrreducibleSolvableGroups( t[1], t[2] )));
372
gap> STOP_TEST( "irrednumbers.tst", 1);
