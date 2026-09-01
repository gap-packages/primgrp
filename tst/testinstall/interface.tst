#@local grps
gap> START_TEST( "interface.tst" );

#
gap> grps:= AllPrimitiveGroups( NrMovedPoints, 7 );;
gap> Length( grps );
7
gap> AllPrimitiveGroups( NrMovedPoints, 7, NrMovedPoints, 8 );
[  ]
gap> AllPrimitiveGroups( NrMovedPoints, IsPrimeInt, NrMovedPoints, 8 );
[  ]
gap> AllPrimitiveGroups( NrMovedPoints, IsPrimeInt, NrMovedPoints, [ 8, 9 ] );
[  ]
gap> AllPrimitiveGroups( NrMovedPoints, IsPrimeInt, NrMovedPoints, [ 6 .. 10 ] ) =
>    grps;
true

#
gap> List( AllPrimitiveGroups( Size, [ 1 .. 100 ], Size, IsPrimeInt ), Size ) =
>    Filtered( [ 1 .. 100 ], IsPrimeInt );
true

#
gap> STOP_TEST( "interface.tst" );
