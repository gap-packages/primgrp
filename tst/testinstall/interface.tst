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

# Only the intersection of the degree lists decides whether the search
# was complete; here it is [ 7 ] resp. [ ], although each single list
# reaches beyond the library.
gap> AllPrimitiveGroups( NrMovedPoints, [ 7, 9000 ],
>                        NrMovedPoints, [ 7, 10000 ] ) = grps;
true
gap> AllPrimitiveGroups( NrMovedPoints, [ 1, 3 .. 8193 ],
>                        NrMovedPoints, [ 2, 4 .. 8192 ] );
[  ]
gap> AllPrimitiveGroups( NrMovedPoints, [ 7, 9000 ] ) = grps;
#W  AllPrimitiveGroups: Degree restricted to [ 2 .. 8191 ]
true

#
gap> List( AllPrimitiveGroups( Size, [ 1 .. 100 ], Size, IsPrimeInt ), Size ) =
>    Filtered( [ 1 .. 100 ], IsPrimeInt );
true

# All `Size` and `Order` conditions restrict the degree, also next to
# an explicit `NrMovedPoints` condition.
gap> AllPrimitiveGroups( Size, 168, Size, 336 );
[  ]
gap> ForAll( AllPrimitiveGroups( Size, [ 1 .. 100 ], Order, [ 50 .. 200 ] ),
>            g -> Size( g ) in [ 50 .. 100 ] );
true
gap> List( AllPrimitiveGroups( NrMovedPoints, [ 1 .. 200 ], Size, 168 ),
>          NrMovedPoints ) =
>    List( AllPrimitiveGroups( Size, 168 ), NrMovedPoints );
true

#
gap> STOP_TEST( "interface.tst" );
