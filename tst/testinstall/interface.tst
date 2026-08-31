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
gap> Length( AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, 0 ) ) = 0;
true

#
gap> STOP_TEST( "interface.tst" );
