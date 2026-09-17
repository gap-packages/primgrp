#############################################################################
##
##  dev/convert-extended.g
##
##  Describe the entries of degree 4096 to 8191 that dev/import-extended.g
##  left as generators.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" \
##          -c 'conv_dir:="data";; conv_first:=40;; conv_4c:="<file>";;' \
##          dev/convert-extended.g
##
##  The product action entries take their wreath elements from a file holding
##  the entries dev/product4c.g wrote for this data, each proven there by
##  rebuilding it and checking it equals the original after relabelling the
##  points.  Field 9 is taken only if fields 1 to 8 are identical, and the
##  group rebuilt from it must have the entry's order, transitivity and
##  suborbits and be primitive.
##
##  The actions of the L series on points, pairs and k-spaces, and of Alt and
##  Sym on k-sets, become descriptions as below degree 4096, but with a proof
##  that avoids PrimitiveIdentification.  The group built must be primitive of
##  the degree with the entry's order, transitivity and suborbits, and no other
##  entry of the degree may have the same order, O'Nan-Scott type, socle,
##  suborbits and transitivity.  Every primitive group of the degree is
##  conjugate to exactly one entry, so that makes it this one.  Where
##  2-transitive entries share all of these, the order and cyclicity of a
##  two-point stabiliser must separate them.  Entries not separated stay as
##  generators and are printed.
##
##  Only invariants computed from the group built enter the fingerprint.  Its
##  order is proven before a stabiliser chain in the large degree is asked for:
##  an action on k-sets or k-spaces is faithful, so its order is the inner
##  group's; a product action has the order of the same generators acting on
##  m*k points.
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

##  Where two of the four coincide the first is taken.
PRIMGRP_Names := ["PSL", "PGL", "PGammaL", "PSigmaL"];
PRIMGRP_Ctor := rec(PSL := PGPsl, PGL := PGPgl, PGammaL := PGPgammaL,
                    PSigmaL := PGPsigmaL);
PRIMGRP_GapGroup := rec(PSL := PSL, PGL := PGL, PGammaL := PGammaL,
                        PSigmaL := PSigmaL);

PRIMGRP_Log := rec(fourc := 0, points := 0, words := 0, sets := 0,
                   subspaces := 0, fourcbad := [], refused := [], pending := [],
                   renamed := [], report := []);

##  A constructor's entry agrees with the stored one in all but name and field 9.
PRIMGRP_Agrees := function(c, e)
  return c{[1..4]} = e{[1..4]} and Set(c[5]) = Set(e[5]) and c[6] = e[6]
         and c[8] = e[8];
end;

PRIMGRP_Fingerprint := e -> [e[2], Set(e[5]), e[6]];

PRIMGRP_TwoPoint := function(g)
  local s;
  s := Stabilizer(Stabilizer(g, 1), 2);
  return [Size(s), IsCyclic(s)];
end;

##  Is g, built for entry nr of degree deg, that entry?
PRIMGRP_Proves := function(g, deg, nr)
  local e, fp, sib, inv;
  e := PRIMGrp(deg, nr);
  if NrMovedPoints(g) <> deg or Size(g) <> e[2]
     or not IsPrimitive(g, [1..deg]) or Transitivity(g, [1..deg]) <> e[6]
     or Set(Collected(OrbitLengthsDomain(Stabilizer(g, 1), [2..deg])))
        <> Set(e[5]) then
    return false;
  fi;
  fp := PRIMGRP_Fingerprint(e);
  sib := Filtered([1..NrPrimitiveGroups(deg)],
                  i -> i <> nr and PRIMGRP_Fingerprint(PRIMGrp(deg, i)) = fp);
  if sib = [] then
    return true;
  fi;
  if e[6] >= 2 then
    inv := PRIMGRP_TwoPoint(g);
    if PRIMGRP_TwoPoint(PrimitiveGroup(deg, nr)) = inv
       and ForAll(sib, i -> PRIMGRP_TwoPoint(PrimitiveGroup(deg, i)) <> inv) then
      return true;
    fi;
  fi;
  AddSet(PRIMGRP_Log.pending, [deg, nr, sib]);
  return false;
end;

PRIMGRP_Renamed := function(deg, e, c)
  if c[7] <> e[7] then
    Add(PRIMGRP_Log.renamed, [deg, e[1], e[7], c[7]]);
  fi;
end;

##  The lists of words generating the subgroups of PGammaL/PSL other than the
##  four named ones, the first list for each subgroup.
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

PRIMGRP_DescribePoints := function(e, deg, dim, q)
  local name, c, words;
  for name in PRIMGRP_Names do
    c := PRIMGRP_Ctor.(name)(dim, q)(deg, e[1]);
    if PRIMGRP_Agrees(c, e)
       and PRIMGRP_Proves(PRIMGRP_GapGroup.(name)(dim, q), deg, e[1]) then
      PRIMGRP_Renamed(deg, e, c);
      PRIMGRP_Log.points := PRIMGRP_Log.points + 1;
      return PRIMGRP_Compact([name, dim, q]);
    fi;
  od;
  for words in PRIMGRP_ExtensionWords(dim, q) do
    c := PGPslExtended(dim, q, words, e[7])(deg, e[1]);
    if PRIMGRP_Agrees(c, e)
       and PRIMGRP_Proves(PGOnPointsGroup(["PSL", dim, q, words, e[7]]),
                          deg, e[1]) then
      PRIMGRP_Log.words := PRIMGRP_Log.words + 1;
      return PRIMGRP_Compact(["PSL", dim, q, words, e[7]]);
    fi;
  od;
  return fail;
end;

##  PSL(2,q) and its companions on pairs of points, and Alt and Sym on k-sets.
PRIMGRP_DescribeSets := function(e, deg, inners, k)
  local inner, c, g;
  for inner in inners do
    c := PGOnSets(inner, k)(deg, e[1]);
    if not PRIMGRP_Agrees(c, e) then
      continue;
    fi;
    g := PGOnSetsGroup(inner, k);
    SetSize(g, Size(PRIMGRP_InnerGroup(inner)));
    if PRIMGRP_Proves(g, deg, e[1]) then
      PRIMGRP_Renamed(deg, e, c);
      PRIMGRP_Log.sets := PRIMGRP_Log.sets + 1;
      return PRIMGRP_Compact(["sets", inner, k]);
    fi;
  od;
  return fail;
end;

PRIMGRP_DescribeSubspaces := function(e, deg, dim, q, k)
  local name, inner, c, g;
  for name in PRIMGRP_Names do
    inner := [name, dim, q];
    c := PGOnSubspaces(inner, k)(deg, e[1]);
    if not PRIMGRP_Agrees(c, e) then
      continue;
    fi;
    g := PGOnSubspacesGroup(inner, k);
    SetSize(g, Size(PRIMGRP_InnerGroup(inner)));
    if PRIMGRP_Proves(g, deg, e[1]) then
      PRIMGRP_Renamed(deg, e, c);
      PRIMGRP_Log.subspaces := PRIMGRP_Log.subspaces + 1;
      return PRIMGRP_Compact(["subspaces", inner, k]);
    fi;
  od;
  return fail;
end;

PRIMGRP_Transplant := rec();

##  The order of the group PGProductAction4c(m, k, els) builds, from its
##  action on k copies of [1..m]: coordinate i holding c goes to coordinate
##  i^sigma holding c^p_i.  The action on tuples is induced from it and both
##  are faithful for m >= 2.
PRIMGRP_WreathOrder := function(m, k, els)
  return Size(Group(List(els, e -> PermList(Concatenation(List([1..k],
           i -> List([1..m], c -> (i^e[k+1] - 1)*m + c^e[i]))))), ()));
end;

PRIMGRP_Describe4c := function(e, deg)
  local key, old, k, H;
  key := Concatenation(String(deg), "_", String(e[1]));
  if not IsBound(PRIMGRP_Transplant.(key)) then
    return fail;
  fi;
  old := PRIMGRP_Transplant.(key);
  if not (PRIMGRP_Agrees(old, e) and old[7] = e[7]) then
    Add(PRIMGRP_Log.fourcbad, [deg, e[1], "fields differ"]);
    return fail;
  fi;
  k := e[8][3];
  H := PGProductAction4c(RootInt(deg, k), k, old[9][2]);
  if RootInt(deg, k) < 2
     or PRIMGRP_WreathOrder(RootInt(deg, k), k, old[9][2]) <> e[2] then
    Add(PRIMGRP_Log.fourcbad, [deg, e[1], "order differs"]);
    return fail;
  fi;
  SetSize(H, e[2]);
  if NrMovedPoints(H) <> deg
     or Transitivity(H, [1..deg]) <> e[6]
     or Set(Collected(OrbitLengthsDomain(Stabilizer(H, 1), [2..deg])))
        <> Set(e[5])
     or not IsPrimitive(H, [1..deg]) then
    Add(PRIMGRP_Log.fourcbad, [deg, e[1], "rebuilt group differs"]);
    return fail;
  fi;
  e := ShallowCopy(e);
  e[9] := old[9];
  PRIMGRP_Log.fourc := PRIMGRP_Log.fourc + 1;
  return PRIMGRP_Compact(e);
end;

##  The text for an entry still holding generators, or fail.
PRIMGRP_Describe := function(e, deg)
  local dim, q, n, k, m, t;
  if e[4] = "4c" then
    return PRIMGRP_Describe4c(e, deg);
  fi;
  if e[4] <> "2" or not (IsList(e[8]) and Length(e[8]) = 3 and e[8][3] = 1) then
    return fail;
  fi;
  if e[8][1] = "L" and IsList(e[8][2]) then
    dim := e[8][2][1];
    q := e[8][2][2];
    n := (q^dim-1)/(q-1);
    if deg = n then
      return PRIMGRP_DescribePoints(e, deg, dim, q);
    fi;
    if dim = 2 and deg = Binomial(n, 2) then
      t := PRIMGRP_DescribeSets(e, deg,
             List(PRIMGRP_Names, name -> [name, dim, q]), 2);
      if t <> fail then
        return t;
      fi;
    fi;
    k := First([2..QuoInt(dim, 2)],
               i -> PRIMGRP_GaussianBinomial(dim, i, q) = deg);
    if k <> fail then
      return PRIMGRP_DescribeSubspaces(e, deg, dim, q, k);
    fi;
  elif e[8][1] = "A" and IsInt(e[8][2]) then
    m := e[8][2];
    k := First([2..QuoInt(m, 2)], i -> Binomial(m, i) = deg);
    if k <> fail then
      return PRIMGRP_DescribeSets(e, deg, [["Alt", m], ["Sym", m]], k);
    fi;
  fi;
  return fail;
end;

PRIMGRP_ConvertFile := function(path)
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
    elif line{[1..2]} = "];" or line{[1..2]} = "[\""
         or line[Length(line)] <> ','
         or (PositionSublist(line, ",\"2\",") = fail
             and PositionSublist(line, ",\"4c\",") = fail) then
      Add(out, line);
    else
      e := EvalString(line{[1..Length(line)-1]});
      if not ForAll(e[9], IsPerm) then
        Add(out, line);
        continue;
      fi;
      # a candidate that fails the proof is a refusal only if none succeeds
      PRIMGRP_Log.pending := [];
      text := PRIMGRP_Describe(e, deg);
      if text = fail then
        UniteSet(PRIMGRP_Log.refused, PRIMGRP_Log.pending);
        Add(out, line);
      else
        Add(out, Concatenation(text, ","));
        Add(PRIMGRP_Log.report, [deg, e[1], text{[1..Minimum(60, Length(text))]}]);
        n := n + 1;
      fi;
    fi;
  od;
  if n > 0 then
    FileString(path, Concatenation(JoinStringsWithSeparator(out, "\n"), "\n"));
  fi;
  return n;
end;

PRIMGRP_ConvertAll := function()
  local line, parts, f, n, i, r;
  for line in SplitString(StringFile(conv_4c), "\n") do
    if Length(line) = 0 then
      continue;
    fi;
    parts := SplitString(line, "\t");
    PRIMGRP_Transplant.(Concatenation(parts[1], "_", parts[2])) :=
      EvalString(parts[3]);
  od;
  n := 0;
  i := conv_first;
  f := Concatenation(conv_dir, "/gps", String(i), ".g");
  while IsReadableFile(f) do
    n := n + PRIMGRP_ConvertFile(f);
    # the proofs load every entry of a degree; free them with the file
    for r in Filtered([4096..8191], d -> IsBound(PRIMGRP[d])) do
      Unbind(PRIMGRP[r]);
    od;
    Print("FILE gps", i, ": ", n, " converted so far\n");
    i := i + 1;
    f := Concatenation(conv_dir, "/gps", String(i), ".g");
  od;
  for r in PRIMGRP_Log.report do
    Print("  ", r[1], "/", r[2], "  ", r[3], "\n");
  od;
  for r in PRIMGRP_Log.renamed do
    Print("RENAMED ", r[1], "/", r[2], "  ", ViewString(r[3]), " -> ",
          ViewString(r[4]), "\n");
  od;
  for r in PRIMGRP_Log.refused do
    Print("REFUSED ", r[1], "/", r[2], "  shares its fingerprint with ",
          r[3], "\n");
  od;
  for r in PRIMGRP_Log.fourcbad do
    Print("4C NOT TAKEN ", r[1], "/", r[2], "  ", r[3], "\n");
  od;
  Print("CONVERTED ", n, " entries: 4c ", PRIMGRP_Log.fourc, ", points ",
        PRIMGRP_Log.points, ", extended ", PRIMGRP_Log.words, ", sets ",
        PRIMGRP_Log.sets, ", subspaces ", PRIMGRP_Log.subspaces, "; renamed ",
        Length(PRIMGRP_Log.renamed), ", refused ", Length(PRIMGRP_Log.refused),
        ", 4c not taken ", Length(PRIMGRP_Log.fourcbad), "\n");
end;

PRIMGRP_ConvertAll();
QUIT;
