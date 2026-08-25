#############################################################################
##
##  dev/find-alt-sym.g
##
##  Find the entries that PGAlt or PGSym reproduces exactly.
##
##  Not "the entry is named Alt(n)" -- the name is the least trustworthy field
##  in the file, and above degree 4095 there is often none at all.  The test is
##  that the constructor's output *equals the entry*, all nine fields, as
##  values rather than as text.  An entry that passes it can be replaced by the
##  call without asking what it was, because the call is already known to
##  produce it.
##
##  That also disposes of the small degrees without a special case: Alt(4) is
##  affine and solvable, so PGAlt(4,1) does not equal it and it is left alone.
##
##  Parameters, via -c:  fs_from, fs_to, fs_out
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);

PRIMGRP_FindAltSym := function(from, to, file)
  local out, deg, nr, l, n, seen;
  out := OutputTextFile(file, false);
  SetPrintFormattingStatus(out, false);
  n := 0;
  seen := 0;
  for deg in [from..to] do
    for nr in [1..NrPrimitiveGroups(deg)] do
      l := PRIMGrp(deg, nr);
      seen := seen + 1;
      if IsList(l) then
        if l = PGAlt(deg, nr) then
          PrintTo(out, deg, " ", nr, " PGAlt\n");
          n := n + 1;
        elif l = PGSym(deg, nr) then
          PrintTo(out, deg, " ", nr, " PGSym\n");
          n := n + 1;
        fi;
      fi;
    od;
  od;
  CloseStream(out);
  Print("FOUND ", n, " of ", seen, " entries reproduced exactly, degrees ",
        from, " to ", to, "\n");
end;

PRIMGRP_FindAltSym(fs_from, fs_to, fs_out);
QUIT;
