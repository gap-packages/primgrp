#@local d
gap> START_TEST( "interface.tst" );

#
gap> for d in [ 1 .. 100 ] do
>      if AllPrimitiveGroups( NrMovedPoints, d, IsAlmostSimple, false ) <>
>      AllPrimitiveGroups( NrMovedPoints, d, ONanScottType, x -> x <> "2" ) then
>        Error( "inconsistency" );
>      fi;
>    od;

#
gap> STOP_TEST( "interface.tst" );
