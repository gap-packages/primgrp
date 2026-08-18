#############################################################################
##
##  dev/regen-bundled.g
##
##  Rewrite every data/gps*.g in the compact format, keeping the existing
##  degree-to-file map so PRIMINDX is untouched.  Writes to PRIMGRP_OutDir;
##  the caller compares before swapping anything in.
##
Read("dev/generic.g");
Read("dev/serialize.g");

BindGlobal("PRIMGRP_RegenBundled", function(outDir)
  local ind, degs, deg, out, name, total;
  total := 0;
  for ind in Set(PRIMINDX{[2..8191]}) do
    degs := Filtered([2..8191], d -> PRIMINDX[d] = ind);
    out := "";
    for deg in degs do
      PrimGrpLoad(deg);
      Append(out, PRIMGRP_CompactDegree(
        List([1..NrPrimitiveGroups(deg)], nr -> PRIMGrp(deg, nr)), deg));
    od;
    name := Concatenation(outDir, "/gps", String(ind), ".g");
    FileString(name, out);
    total := total + Length(out);
  od;
  return total;
end);
