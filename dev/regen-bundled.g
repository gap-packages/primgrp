#############################################################################
##
##  dev/regen-bundled.g
##
##  Rewrite data/gps1.g .. gps39.g in the compact format, keeping the existing
##  degree-to-file map so PRIMINDX is untouched.  Writes to PRIMGRP_OutDir;
##  the caller compares before swapping anything in.
##
Read("dev/generic.g");
Read("dev/serialize.g");

BindGlobal("PRIMGRP_RegenBundled", function(outDir)
  local ind, degs, deg, out, name, total;
  total := 0;
  for ind in Set(PRIMINDX{[2..4095]}) do
    degs := Filtered([2..4095], d -> PRIMINDX[d] = ind);
    out := "";
    for deg in degs do
      PrimGrpLoad(deg);
      Append(out, PRIMGRP_CompactDegree(PRIMGRP[deg], deg));
    od;
    name := Concatenation(outDir, "/gps", String(ind), ".g");
    FileString(name, out);
    total := total + Length(out);
  od;
  return total;
end);
