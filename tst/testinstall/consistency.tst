gap> START_TEST("consistency.tst");

#
# Warnings here depend on the Conway polynomial database
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
gap> ReadPackage("primgrp", "tst/testutils.g");
true

#
# Size, Transitivity, ONanScottType, IsSimpleGroup, IsSolvableGroup, the
# collected suborbits and SocleTypePrimitiveGroup are all served from the data
# files rather than computed.  Recompute them and compare.  The full sweep over
# every degree lives in tst/testextra/consistency.tst.
#
gap> Concatenation(List([2..100], d -> PrimGrpCheckDegree(d, PrimGrpCheckAll)));
[  ]

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "consistency.tst", 1);
