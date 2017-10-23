gap> START_TEST("primid.tst");
gap> ForAll([2..4095], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> STOP_TEST( "primid.tst", 1);
