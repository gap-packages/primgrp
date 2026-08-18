#############################################################################
##
##  dev/check-sigma-gamma.g
##
##  Check PSigmaL(dim,q) and PGammaL(dim,q) against the entries they would
##  replace: degree, order, transitivity, suborbits, primitivity.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_CheckSG", function(items, logFile)
  local it, kind, deg, nr, dim, q, e, H, t, out, bad, line;
  for it in items do
    kind := it[1]; deg := it[2]; nr := it[3]; dim := it[4]; q := it[5];
    e := PRIMGrp(deg, nr);
    if kind = "sigma" then H := PSigmaL(dim, q); else H := PGammaL(dim, q); fi;
    bad := [];
    if NrMovedPoints(H) <> deg then Add(bad, "degree"); fi;
    if Size(H) <> e[2] then
      Add(bad, Concatenation("order stored ", String(e[2]),
                             " computed ", String(Size(H))));
    fi;
    t := Transitivity(H, [1..deg]);
    if t <> e[6] then
      Add(bad, Concatenation("transitivity stored ", String(e[6]),
                             " computed ", String(t)));
    fi;
    if PrimGrpSuborbits(H, deg) <> Set(e[5]) then Add(bad, "suborbits"); fi;
    if not IsPrimitive(H, [1..deg]) then Add(bad, "not primitive"); fi;
    line := Concatenation(kind, " ", String(deg), " ", String(nr),
                          " (", String(dim), ",", String(q), ") ");
    if IsEmpty(bad) then Append(line, "ok\n");
    else Append(line, Concatenation("BAD ", JoinStringsWithSeparator(bad, "; "), "\n")); fi;
    out := OutputTextFile(logFile, true);
    SetPrintFormattingStatus(out, false);
    WriteAll(out, line);
    CloseStream(out);
  od;
end);

PRIMGRP_CheckSG(PRIMGRP_Items, PRIMGRP_LogFile);
QUIT;
