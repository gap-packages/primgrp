gap> START_TEST("bugfix.tst");

# Use consistent names for the Tits group and for the Ree group containing it.
# See <https://github.com/gap-packages/primgrp/issues/17>
gap> G := PrimitiveGroup(1600, 20);
2F(4, 2)'
gap> IsomorphismTypeInfoFiniteSimpleGroup(G);
rec( name := "2F(4,2)' = Ree(2)' = Tits", parameter := 2, series := "2F", 
  shortname := "2F4(2)'" )

#
gap> G := PrimitiveGroup(1755, 1);
2F(4, 2)'
gap> IsomorphismTypeInfoFiniteSimpleGroup(G);
rec( name := "2F(4,2)' = Ree(2)' = Tits", parameter := 2, series := "2F", 
  shortname := "2F4(2)'" )

#
gap> G := PrimitiveGroup(1755, 2);
2F(4, 2)

#
gap> G := PrimitiveGroup(2304, 9);
2F(4, 2)'
gap> IsomorphismTypeInfoFiniteSimpleGroup(G);
rec( name := "2F(4,2)' = Ree(2)' = Tits", parameter := 2, series := "2F", 
  shortname := "2F4(2)'" )

#
gap> G := PrimitiveGroup(2304, 10);
2F(4, 2)

#
gap> STOP_TEST("bugfix.tst", 1);
