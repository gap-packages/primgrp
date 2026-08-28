#@local d
gap> START_TEST( "interface.tst" );

#
gap> for d in [ 1 .. 100 ] do
>      if AllPrimitiveGroups( NrMovedPoints, d, IsAlmostSimple, false ) <>
>      AllPrimitiveGroups( NrMovedPoints, d, ONanScottType, x -> x <> "2" ) then
>        Error( "inconsistency" );
>      fi;
>    od;
gap> AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, true ) =
>    AllPrimitiveGroups( NrMovedPoints, 8, IsAlmostSimple, x -> x <> false );
true

#
gap> STOP_TEST( "interface.tst" );
