#############################################################################
##
##  dev/product-action.g
##
##  Type 4c is the product action: the socle is T^k and the group sits inside
##  H wr P acting on Delta^k, where H is primitive of degree m on Delta and P
##  is transitive of degree k.  The stored entry already gives k, as the width
##  of the socle type, and m, because the degree is exactly m^k.
##
##  So rather than searching per entry, enumerate the candidate (m,b,k,t) and
##  ask PrimitiveIdentification which library entry each one is.  There are
##  only 1, 2, 5 and 5 transitive groups of degree 2, 3, 4 and 5, so the space
##  is small; what is expensive is the identification itself, hence the
##  per-combination marker files.
##
LoadPackage("primgrp");
LoadPackage("transgrp");
SetInfoLevel(InfoWarning, 0);

if not IsBound(PRIMGRP_MaxDegree) then PRIMGRP_MaxDegree := 8191; fi;

BindGlobal("PRIMGRP_ScanProductAction", function(bases, markDir, logFile)
  local m, b, k, t, deg, W, id, mark, out, line;
  for m in bases do
    for k in [2..5] do
      deg := m^k;
      if deg > PRIMGRP_MaxDegree then continue; fi;
      for b in [1..NrPrimitiveGroups(m)] do
        for t in [1..NrTransitiveGroups(k)] do
          mark := Concatenation(markDir, "/w", String(m), "_", String(b), "_",
                                String(k), "_", String(t));
          if IsReadableFile(mark) then continue; fi;
          FileString(mark, "");
          W := WreathProductProductAction(PrimitiveGroup(m,b), TransitiveGroup(k,t));
          line := Concatenation(String(m), " ", String(b), " ", String(k), " ",
                                String(t), " ", String(deg), " ");
          if not IsPrimitive(W, [1..deg]) then
            Append(line, "notprimitive\n");
          else
            id := PrimitiveIdentification(W);
            Append(line, Concatenation("is ", String(id), " order ",
                                       String(Size(W)), "\n"));
          fi;
          out := OutputTextFile(logFile, true);
          SetPrintFormattingStatus(out, false);
          WriteAll(out, line);
          CloseStream(out);
        od;
      od;
    od;
  od;
end);

PRIMGRP_ScanProductAction(PRIMGRP_Bases, PRIMGRP_MarkDir, PRIMGRP_LogFile);
QUIT;
