#############################################################################
##
##  dev/convert-points-projective.g
##
##  Write the actions of the L series on points as descriptions.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'conv_dir:="data";;' \
##          dev/convert-points-projective.g
##
##  An entry whose socle is L(dim,q) and whose degree is the number of points
##  of PG(dim-1,q) lies between PSL(dim,q) and PGammaL(dim,q).
##
##  If it is one of PSL, PGL, PSigmaL and PGammaL it becomes that description
##  and takes the constructor's name; renames are printed.  Where two coincide
##  the first is taken, in the order PSL, PGL, PGammaL, PSigmaL.  Otherwise it
##  becomes ["PSL",dim,q,<words>,<name>], keeping its name, for the first list
##  of words generating a subgroup of PGammaL/PSL of the entry's index.
##
##  Either way the constructor must reproduce every other field, and
##  PrimitiveIdentification must name the group the reader builds as that
##  entry.
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

##  Name, entry constructor, and the group the reader builds for that entry.
PRIMGRP_Named := [ ["PSL", PGPsl, PSL], ["PGL", PGPgl, PGL],
                   ["PGammaL", PGPgammaL, PGammaL],
                   ["PSigmaL", PGPsigmaL, PSigmaL] ];

##  The lists of words generating the subgroups of PGammaL/PSL other than
##  the four named ones, the first list for each subgroup.
PRIMGRP_ExtensionWords := function(dim, q)
  local f, g, els, named, seen, out, gens, H;
  f := Length(Factors(q));
  g := Gcd(dim, q-1);
  els := Difference(Cartesian([0..g-1], [0..f-1]), [[0,0]]);
  named := List([[], [[1 mod g,0]], [[0,1 mod f]], [[1 mod g,0],[0,1 mod f]]],
                w -> PRIMGRP_ProjectiveExtension(dim, q, w));
  seen := ShallowCopy(named);
  out := [];
  for gens in Concatenation(List(els, x -> [x]), Combinations(els, 2)) do
    H := PRIMGRP_ProjectiveExtension(dim, q, gens);
    if not H in seen then
      Add(seen, H);
      Add(out, gens);
    fi;
  od;
  return out;
end;

PRIMGRP_IdCache := rec();

PRIMGRP_Identify := function(key, build, deg)
  local g;
  if not IsBound(PRIMGRP_IdCache.(key)) then
    g := build();
    if NrMovedPoints(g) <> deg then
      PRIMGRP_IdCache.(key) := fail;
    else
      PRIMGRP_IdCache.(key) := PrimitiveIdentification(g);
    fi;
  fi;
  return PRIMGRP_IdCache.(key);
end;

PRIMGRP_Agrees := function(c, e)
  return c{[1..4]} = e{[1..4]} and Set(c[5]) = Set(e[5]) and c[6] = e[6]
         and c[8] = e[8];
end;

PRIMGRP_ConvertPoints := function(e, deg, renamed)
  local dim, q, i, c, key, words;
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
  if deg <> (q^dim-1)/(q-1) then
    return fail;
  fi;

  for i in [1..4] do
    c := PRIMGRP_Named[i][2](dim, q)(deg, e[1]);
    key := Concatenation(PRIMGRP_Named[i][1], "_", String(dim), "_", String(q));
    if PRIMGRP_Agrees(c, e)
       and PRIMGRP_Identify(key, {} -> PRIMGRP_Named[i][3](dim, q), deg)
           = e[1] then
      if c[7] <> e[7] then
        Add(renamed, [deg, e[1], e[7], c[7]]);
      fi;
      return PRIMGRP_Compact([PRIMGRP_Named[i][1], dim, q]);
    fi;
  od;

  for words in PRIMGRP_ExtensionWords(dim, q) do
    c := PGPslExtended(dim, q, words, e[7])(deg, e[1]);
    key := Concatenation("PSL_", String(dim), "_", String(q), "_",
             JoinStringsWithSeparator(List(words, w ->
               Concatenation(String(w[1]), "x", String(w[2]))), "_"));
    if PRIMGRP_Agrees(c, e)
       and PRIMGRP_Identify(key,
             {} -> PGOnPointsGroup(["PSL", dim, q, words]), deg) = e[1] then
      return PRIMGRP_Compact(["PSL", dim, q, words, e[7]]);
    fi;
  od;
  return fail;
end;

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
    elif PositionSublist(line, ",[\"L\",[") = fail then
      Add(out, line);                   # cheap test before EvalString
    else
      e := EvalString(line{[1..Length(line)-1]});
      text := PRIMGRP_ConvertPoints(e, deg, renamed);
      if text = fail then
        Add(out, line);
      else
        Add(out, Concatenation(text, ","));
        Add(report, [deg, e[1], text]);
        n := n + 1;
      fi;
    fi;
  od;
  if n > 0 then
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
  Print("CONVERTED ", n, " entries, ", Length(renamed), " renamed\n");
end;

PRIMGRP_ConvertAll(conv_dir);
QUIT;
