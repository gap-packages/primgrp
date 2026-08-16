gap> START_TEST("consistency.tst");

#
# The whole library, stored against recomputed.  This takes hours: it builds a
# stabiliser chain for every group that is not a natural A_n or S_n, and there
# are 54334 groups.  Not part of any automatic test run.
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
gap> ReadPackage("primgrp", "tst/testutils.g");
true

#
gap> Concatenation(List([2..8191], d -> PrimGrpCheckDegree(d, PrimGrpCheckAll)));
[  ]

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "consistency.tst", 1);
