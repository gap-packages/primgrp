#@local d, grps
gap> START_TEST( "interface.tst" );

#
gap> for d in [ 1 .. 100 ] do
>      if AllPrimitiveGroups( NrMovedPoints, d, IsAlmostSimple, false ) <>
>      AllPrimitiveGroups( NrMovedPoints, d, ONanScottType, x -> x <> "2" ) then
>        Error( "inconsistency" );
>      fi;
>    od;
gap> grps:= AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, true );;
gap> grps =
>    AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, x -> x <> false );
true
gap> grps =
>    AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, [ true ] );
true
gap> AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, [ true, false ] ) =
>    AllPrimitiveGroups( NrMovedPoints, 8 );
true
gap> Length( AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, fail ) ) = 0;
true

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
