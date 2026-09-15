#############################################################################
##
##  dev/convert-sets-projective.g
##
##  Write the actions on k-subsets as ["sets", <inner>, k].
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'conv_dir:="data";;' \
##          dev/convert-sets-projective.g
##
##  ["AltOnSets",n,k] and ["SymOnSets",n,k] are respelled as
##  ["sets", ["Alt",n], k] and ["sets", ["Sym",n], k], which build the same
##  entry.
##
##  An entry whose socle is L(dim,q) and whose degree is Binomial(n,k), n the
##  number of points of PG(dim-1,q), may be PSL, PGL, PSigmaL or PGammaL on
##  the k-subsets of the points.  It then becomes ["sets", ["PSL",dim,q], k]
##  and so on, provided the constructor reproduces every field but the name and
##  PrimitiveIdentification names the group the reader builds as that entry.
##  The name becomes the constructor's; renames are printed.  Where two of the
##  four coincide the first is taken, in the order PSL, PGL, PGammaL, PSigmaL.
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

PRIMGRP_Compact := function(o)
  if IsStringRep(o) then
    return ViewString(o);
  fi;
  if IsList(o) then
    return Concatenation("[", JoinStringsWithSeparator(
             List(o, PRIMGRP_Compact), ","), "]");
  fi;
  return String(o);
end;

PRIMGRP_Names := ["PSL", "PGL", "PGammaL", "PSigmaL"];

PRIMGRP_IdCache := rec();

PRIMGRP_Identify := function(inner, k, deg)
  local key, g;
  key := Concatenation(inner[1], "_", String(inner[2]), "_", String(inner[3]),
                       "_", String(k));
  if not IsBound(PRIMGRP_IdCache.(key)) then
    g := PGOnSetsGroup(inner, k);
    if NrMovedPoints(g) <> deg then
      PRIMGRP_IdCache.(key) := fail;
    else
      PRIMGRP_IdCache.(key) := PrimitiveIdentification(g);
    fi;
  fi;
  return PRIMGRP_IdCache.(key);
end;

PRIMGRP_ConvertSets := function(e, deg, renamed)
  local dim, q, n, k, name, c;
  if not (IsList(e) and Length(e) = 9 and e[4] = "2") then
    return fail;
  fi;
  if not (IsList(e[8]) and Length(e[8]) = 3 and e[8][1] = "L"
          and IsList(e[8][2]) and e[8][3] = 1) then
    return fail;
  fi;
  if IsString(e[9]) or ForAny(e[9], x -> not IsPerm(x)) then
    return fail;
  fi;
  dim := e[8][2][1];
  q := e[8][2][2];
  n := (q^dim-1)/(q-1);
  k := First([2..QuoInt(n, 2)], i -> Binomial(n, i) = deg);
  if k = fail then
    return fail;
  fi;
  for name in PRIMGRP_Names do
    c := PGOnSets([name, dim, q], k)(deg, e[1]);
    if c{[1..4]} = e{[1..4]} and Set(c[5]) = Set(e[5]) and c[6] = e[6]
       and c[8] = e[8] and PRIMGRP_Identify([name, dim, q], k, deg) = e[1] then
      if c[7] <> e[7] then
        Add(renamed, [deg, e[1], e[7], c[7]]);
      fi;
      return PRIMGRP_Compact(["sets", [name, dim, q], k]);
    fi;
  od;
  return fail;
end;

respelled := 0;

PRIMGRP_ConvertFile := function(path, report, renamed)
  local lines, out, deg, line, e, text, n;
  lines := SplitString(StringFile(path), "\n");
  out := [];
  deg := fail;
  n := 0;
  for line in lines do
    if Length(line) = 0 then
      continue;
    fi;
    if PositionSublist(line, "PRIMGRP[") = 1 then
      deg := Int(line{[9..PositionSublist(line, "]") - 1]});
      Add(out, line);
    elif line{[1..2]} = "];" or line[Length(line)] <> ',' then
      Add(out, line);
    elif PositionSublist(line, "[\"AltOnSets\",") = 1
         or PositionSublist(line, "[\"SymOnSets\",") = 1 then
      e := EvalString(line{[1..Length(line)-1]});
      Add(out, Concatenation(PRIMGRP_Compact(["sets",
                 [e[1]{[1..3]}, e[2]], e[3]]), ","));
      respelled := respelled + 1;
    elif PositionSublist(line, ",[\"L\",[") = fail then
      Add(out, line);                   # cheap test before EvalString
    else
      e := EvalString(line{[1..Length(line)-1]});
      text := PRIMGRP_ConvertSets(e, deg, renamed);
      if text = fail then
        Add(out, line);
      else
        Add(out, Concatenation(text, ","));
        Add(report, [deg, e[1], text]);
        n := n + 1;
      fi;
    fi;
  od;
  if n > 0 or respelled > 0 then
    FileString(path, Concatenation(JoinStringsWithSeparator(out, "\n"), "\n"));
  fi;
  return n;
end;

PRIMGRP_ConvertAll := function(dir)
  local files, f, n, report, renamed, r;
  files := Filtered(SortedList(DirectoryContents(dir)),
                    f -> Length(f) > 5 and f{[1..3]} = "gps"
                         and f{[Length(f)-1..Length(f)]} = ".g");
  n := 0;
  report := [];
  renamed := [];
  for f in files do
    n := n + PRIMGRP_ConvertFile(Concatenation(dir, "/", f), report, renamed);
  od;
  for r in report do
    Print("  ", r[1], "/", r[2], "  ", r[3], "\n");
  od;
  for r in renamed do
    Print("RENAMED ", r[1], "/", r[2], "  ", ViewString(r[3]), " -> ",
          ViewString(r[4]), "\n");
  od;
  Print("RESPELLED ", respelled, " entries\n");
  Print("CONVERTED ", n, " entries, ", Length(renamed), " renamed\n");
end;

PRIMGRP_ConvertAll(conv_dir);
QUIT;
