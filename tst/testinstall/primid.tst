gap> START_TEST("primid.tst");
gap> ForAll([2..1000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([1001..2000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([2001..3000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([3001..4000], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> ForAll([4001..4095], n -> PrimitiveIdentification(OnePrimitiveGroup(NrMovedPoints,n))=1);
true
gap> STOP_TEST( "primid.tst", 1);
