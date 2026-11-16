gap> START_TEST("checknames.tst");

#
# Disable warnings which depend on Conway Polynomial databases
#
gap> iW := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);

#
gap> ReadPackage("primgrp", "lib/testutils.g");
true
gap> PrimGrpNamesCheckDegree(1);
gap> PrimGrpNamesCheckDegree(2);
gap> PrimGrpNamesCheckDegree(3);
gap> PrimGrpNamesCheckDegree(4);
gap> PrimGrpNamesCheckDegree(5);
gap> PrimGrpNamesCheckDegree(6);
gap> PrimGrpNamesCheckDegree(7);
gap> PrimGrpNamesCheckDegree(8);

#
gap> SetInfoLevel(InfoWarning,iW);
gap> STOP_TEST( "checknames.tst", 1);
