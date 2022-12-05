LoadPackage("PrimGrp");
TestDirectory( [
  DirectoriesPackageLibrary("PrimGrp","tst/manualexamples"),
  DirectoriesPackageLibrary("PrimGrp","tst/testinstall"),
  DirectoriesPackageLibrary("PrimGrp","tst/teststandard"),
  ], rec(exitGAP := true));
FORCE_QUIT_GAP(1);
