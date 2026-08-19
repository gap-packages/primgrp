#############################################################################
##
##  dev/classical.g
##
##  Check that PGClassicalGroup rebuilds the type 2 entries listed in
##  dev/classical-unambiguous.txt: those acting on the natural geometry of
##  their socle, whose order is unique among all entries of their degree.
##
##  Uniqueness of the order is the whole argument.  The library is complete
##  up to conjugacy, so a primitive group of that degree with that socle and
##  that order must be conjugate to the one entry that has them -- conjugate,
##  not equal, which is why the check is on invariants rather than on "=".
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_CheckClassical", function(inFile, outFile, logFile)
  local out, log, l, p, deg, nr, ser, d, q, order, ext, g, e, why, ok, n, r;
  out := OutputTextFile(outFile, false);  SetPrintFormattingStatus(out, false);
  log := OutputTextFile(logFile, false);  SetPrintFormattingStatus(log, false);
  n := 0;
  for l in SplitString(StringFile(inFile), "\n") do
    if Length(l) = 0 or l[1] = '#' then continue; fi;
    p := SplitString(l, " ");
    deg := Int(p[1]);  nr := Int(p[2]);  ser := p[3];
    d := Int(p[4]);  q := Int(p[5]);  order := Int(p[6]);  ext := Int(p[8]);
    e := PRIMGrp(deg, nr);
    why := "";  g := fail;
    # Ask what the point action realises before building: an entry whose ext
    # exceeds that is not recoverable this way and must keep its permutations.
    BreakOnError := false;
    r := CALL_WITH_CATCH(PGClassicalNormaliser, [ser, d, q]);
    if not r[1] then
      why := "normaliser failed";
    elif (Size(r[2].A) / Size(r[2].T)) mod ext <> 0 then
      why := Concatenation("outer group is only ",
                           String(Size(r[2].A)/Size(r[2].T)), ", need ", String(ext));
    else
      r := CALL_WITH_CATCH(PGClassicalGroup, [ser, d, q, ext]);
      if r[1] then g := r[2]; else why := "construction failed"; fi;
    fi;
    BreakOnError := true;
    ok := false;
    if g <> fail then
      if   NrMovedPoints(g) <> deg then why := "degree";
      elif Size(g) <> order        then why := "order";
      elif not IsPrimitive(g, [1..deg]) then why := "not primitive";
      elif Transitivity(g, [1..deg]) <> e[6] then why := "transitivity";
      elif PrimGrpSuborbits(g, deg) <> Set(e[5]) then why := "suborbits";
      else ok := true; fi;
    fi;
    n := n + 1;
    if ok then
      PrintTo(out, deg, " ", nr, " [\"cl\",\"", ser, "\",", d, ",", q, ",", ext, "]\n");
      PrintTo(log, deg, " ", nr, " ", ser, "(", d, ",", q, ") ext ", ext, " OK\n");
    else
      PrintTo(log, deg, " ", nr, " ", ser, "(", d, ",", q, ") ext ", ext,
              " SKIP ", why, "\n");
    fi;
  od;
  CloseStream(out);  CloseStream(log);
  Print("CLASSICAL scanned ", n, "\n");
end);
