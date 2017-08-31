LoadPackage("PrimGrp");
TestDirectory( [
  DirectoriesPackageLibrary("PrimGrp","tst/manualexamples"),
  DirectoriesPackageLibrary("PrimGrp","tst/testinstall"),
  ], rec(exitGAP := true));
FORCE_QUIT_GAP(1);
