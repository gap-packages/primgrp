#############################################################################
##
##  dev/check-classical.g
##
##  Check only the socle types with a two-parameter series -- L, B, C, D, 2A,
##  2D, G, 2F -- in degrees >= 2500.
##
##  The general sweep checks everything and so takes many hours; the defect
##  being chased touches only these.  1842 entries rather than 30000 groups.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_ClassicalEntries", function(degrees)
  local out, deg, nr, st;
  out := [];
  for deg in degrees do
    PrimGrpLoad(deg);
    for nr in [1..NrPrimitiveGroups(deg)] do
      st := PRIMGrp(deg, nr)[8];
      if IsList(st[2]) and not st[1] in ["A", "Z"] then
        Add(out, [deg, nr]);
      fi;
    od;
  od;
  return out;
end);

BindGlobal("PRIMGRP_CheckClassical", function(items, markDir, logFile)
  local it, deg, nr, G, st, r, mark, out, line;
  for it in items do
    deg := it[1]; nr := it[2];
    mark := Concatenation(markDir, "/c", String(deg), "_", String(nr));
    if IsReadableFile(mark) then
      continue;
    fi;
    G := PrimitiveGroup(deg, nr);
    st := SocleTypePrimitiveGroup(G);
    r := SocleTypePrimitiveGroup(Group(GeneratorsOfGroup(G)));
    if not PrimGrpSameSocleType(st, r) then
      # Print the width as well.  Leaving it out once cost an afternoon: the
      # socle of PrimitiveGroup(3969,9) has width 2, and a fix written for
      # width 1 matched nothing.
      line := Concatenation("MISMATCH ", String(deg), " ", String(nr),
                            " stored=", st.series, String(st.parameter),
                            " w", String(st.width),
                            " computed=", r.series, String(r.parameter),
                            " w", String(r.width), "\n");
      out := OutputTextFile(logFile, true);
      SetPrintFormattingStatus(out, false);
      WriteAll(out, line);
      CloseStream(out);
    fi;
    FileString(mark, "");
  od;
end);

PRIMGRP_CheckClassical(PRIMGRP_ClassicalEntries(PRIMGRP_Degrees),
                       PRIMGRP_MarkDir, PRIMGRP_LogFile);
QUIT;
