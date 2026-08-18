#############################################################################
##
##  dev/check-pgl.g
##
##  Before rewriting an entry as PGPgl(dim,q), check that PGL(dim,q) really is
##  that group: same order, same transitivity, same suborbits, same socle type.
##
##  Conjugacy again comes from the library being complete, but only for the
##  entries where no other entry of that degree shares this socle and order.
##  That is not automatic here: at degree 170 the entries PGL(2,169),
##  PSigmaL(2,169) and PSL(2,169).2_3 all have the same socle and order and
##  are different groups.  Only the 504 unambiguous ones are passed in.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_CheckPgl", function(items, markDir, logFile)
  local it, deg, nr, dim, q, e, H, t, mark, out, bad, line;
  for it in items do
    deg := it[1]; nr := it[2]; dim := it[3]; q := it[4];
    mark := Concatenation(markDir, "/p", String(deg), "_", String(nr));
    if IsReadableFile(mark) then
      continue;
    fi;
    e := PRIMGrp(deg, nr);
    H := PGL(dim, q);
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

PRIMGRP_CheckPgl(PRIMGRP_Items, PRIMGRP_MarkDir, PRIMGRP_LogFile);
QUIT;
