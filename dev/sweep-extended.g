#############################################################################
##
##  dev/sweep-extended.g
##
##  Convert every degree in 4096-8191, one staging file per degree, and log
##  what happened.  Binning the staging files into data/gps*.g and extending
##  PRIMINDX is a separate, cheap second pass.
##
##  Usage:
##    PRIMGRP_ExtDir := "..../ExtendedPrimitiveGroupsData";;
##    PRIMGRP_StageDir := "..../stage";;
##    PRIMGRP_Degrees := [4096..8191];;
##    Read("dev/sweep-extended.g");
##

Read("dev/convert-extended.g");

BindGlobal("PRIMGRP_Sweep", function(degrees, stageDir, logFile)
  local deg, t, res, s, name, total, tally, line, out;
  total := rec(groups := 0, affine := 0, generic := 0, renamed := 0, perms := 0,
               ffe := 0, int := 0, failed := 0, bytes := 0);
  for deg in degrees do
    name := Concatenation(stageDir, "/deg", String(deg), ".g");
    if IsReadableFile(name) then
      continue;                 # already converted; makes a restart cheap
    fi;
    t := Runtime();
    res := PRIMGRP_ConvertDegree(deg);
    s := PRIMGRP_CompactDegree(res.entries, deg);
    FileString(name, s);
    tally := res.stats;
    total.groups  := total.groups  + Length(res.entries);
    total.affine  := total.affine  + tally.affine;
    total.generic := total.generic + tally.generic;
    total.renamed := total.renamed + tally.renamed;
    total.perms   := total.perms   + tally.perms;
    total.ffe     := total.ffe     + tally.ffe;
    total.int     := total.int     + tally.int;
    total.failed  := total.failed  + Length(tally.failed);
    total.bytes   := total.bytes   + Length(s);
    line := Concatenation(
      String(deg), " groups=", String(Length(res.entries)),
      " affine=", String(tally.affine),
      " ffe=", String(tally.ffe), " int=", String(tally.int),
      " generic=", String(tally.generic),
      " renamed=", String(tally.renamed),
      " perms=", String(tally.perms),
      " failed=", String(Length(tally.failed)),
      " bytes=", String(Length(s)),
      " ms=", String(Runtime()-t), "\n");
    out := OutputTextFile(logFile, true);
    SetPrintFormattingStatus(out, false);
    WriteAll(out, line);
    CloseStream(out);
    if Length(tally.failed) > 0 then
      out := OutputTextFile(logFile, true);
      SetPrintFormattingStatus(out, false);
      WriteAll(out, Concatenation("  FAILED ", String(tally.failed), "\n"));
      CloseStream(out);
    fi;
  od;
  AppendTo(logFile, "TOTAL ", total, "\n");
  return total;
end);

PRIMGRP_Sweep(PRIMGRP_Degrees, PRIMGRP_StageDir, PRIMGRP_LogFile);
