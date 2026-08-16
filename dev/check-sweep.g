#############################################################################
##
##  dev/check-sweep.g
##
##  Run the stored-versus-recomputed check over a range of degrees, logging
##  one line per degree and every complaint.  Resumable: a degree whose marker
##  file exists is skipped.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_CheckSweep", function(degrees, markDir, logFile)
  local deg, t, bad, mark, out, line, x;
  for deg in degrees do
    mark := Concatenation(markDir, "/deg", String(deg));
    if IsReadableFile(mark) then
      continue;
    fi;
    t := Runtime();
    bad := PrimGrpCheckDegree(deg, PrimGrpCheckAll);
    line := Concatenation(String(deg), " groups=", String(NrPrimitiveGroups(deg)),
                          " bad=", String(Length(bad)),
                          " ms=", String(Runtime()-t), "\n");
    for x in bad do
      Append(line, Concatenation("  BAD ", x, "\n"));
    od;
    out := OutputTextFile(logFile, true);
    SetPrintFormattingStatus(out, false);
    WriteAll(out, line);
    CloseStream(out);
    FileString(mark, "");
  od;
end);

PRIMGRP_CheckSweep(PRIMGRP_Degrees, PRIMGRP_MarkDir, PRIMGRP_LogFile);
QUIT;
