gap> START_TEST("cohorts.tst");

#
# Disable warnings which depend on Conway Polynomial databases
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
# Define a function to check the primitive groups of degree n
#
gap> ReadPackage("primgrp", "tst/testutils.g");
true
gap> Perform([2..Length(COHORTS_PRIMITIVE_GROUPS)], n -> CheckCohortConsistency);

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "cohorts.tst", 1);
