#############################################################################
##
##  dev/check-psl.g
##
##  Before rewriting an entry as PGPsl(dim,q), check that PSL(dim,q) really is
##  that group: same order, same transitivity, same suborbits, same socle type.
##
##  Conjugacy comes for free from the library being complete -- no other entry
##  of that degree has this socle and order, and PSL(dim,q) is primitive of
##  that degree, so it can only be conjugate to this one.  That is also what
##  rules out handing back the action on hyperplanes instead of on points.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_CheckPsl", function(items, markDir, logFile)
  local it, deg, nr, dim, q, e, H, t, mark, out, bad, line;
  for it in items do
    deg := it[1]; nr := it[2]; dim := it[3]; q := it[4];
    mark := Concatenation(markDir, "/p", String(deg), "_", String(nr));
    if IsReadableFile(mark) then
      continue;
    fi;
    e := PRIMGrp(deg, nr);
    H := PSL(dim, q);
    bad := [];
    if NrMovedPoints(H) <> deg then Add(bad, "degree"); fi;
    if Size(H) <> e[2] then Add(bad, "size"); fi;
    t := Transitivity(H, [1..deg]);
    if t <> e[6] then
      Add(bad, Concatenation("transitivity stored ", String(e[6]),
                             " computed ", String(t)));
    fi;
    if PrimGrpSuborbits(H, deg) <> Set(e[5]) then Add(bad, "suborbits"); fi;
    if not IsPrimitive(H, [1..deg]) then Add(bad, "not primitive"); fi;
    line := Concatenation(String(deg), " ", String(nr), " L(", String(dim), ",",
                          String(q), ") ");
    if IsEmpty(bad) then
      Append(line, "ok\n");
    else
      Append(line, Concatenation("BAD ", JoinStringsWithSeparator(bad, "; "), "\n"));
    fi;
    out := OutputTextFile(logFile, true);
    SetPrintFormattingStatus(out, false);
    WriteAll(out, line);
    CloseStream(out);
    FileString(mark, "");
  od;
end);

PRIMGRP_CheckPsl(PRIMGRP_Items, PRIMGRP_MarkDir, PRIMGRP_LogFile);
QUIT;
