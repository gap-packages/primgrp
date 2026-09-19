#############################################################################
##
##  dev/import-extended.g
##
##  Bring the primitive groups of degree 4096 to 8191 into data/.
##
##      gap -q -b -A --quitonbreak -l "ROOT;" -c 'imp_archive:="<dir>";;
##          imp_dir:="data";; imp_first:=40;;' \
##          dev/import-extended.g
##
##  It writes the data files and extends PRIMINDX in imp_grp, by default
##  lib/primitiv.grp, to say which file each new degree is in.  Without that
##  the reader refuses the degree, PRIMRANGE or no PRIMRANGE.
##
##  The source is ExtendedPrimitiveGroupsData,
##  https://doi.org/10.5281/zenodo.10411366, one gzipped record per group.  The
##  copy used is the tar with sha256
##  f4379176ae6c271797c2c2dda0799cbf234c1d94b6c79bd1e7fb11ddc6ca1b5c.
##
##  Each entry is written the way the data below 4096 is: the natural
##  alternating and symmetric groups as ["Alt"] and ["Sym"], an affine group of
##  prime degree as ["Prime",d], any other affine group by its matrices, and
##  everything else by its generators.  dev/convert-extended.g then describes
##  what it can of the rest.
##
##  The matrices of an affine group are read off its generators in the
##  archive's numbering of GF(p)^d, where the index of a vector is its base-p
##  digits, and every generator is checked to be the affine map read off it.
##  The reader numbers the points by Elements(GF(p)^d) instead.  The two agree
##  for p = 2 and p = 3; for p >= 5, as for every ["Prime",d], the group read
##  back is a conjugate of the archive's.
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

if not IsBound(imp_perfile) then
  imp_perfile := 800;
fi;
if not IsBound(imp_grp) then
  imp_grp := "lib/primitiv.grp";
fi;

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

##  One record of the archive, as an entry.
PRIMGRP_ReadRecord := function(deg, nr)
  local path, strm, r;
  path := Concatenation(imp_archive, "/PrimitiveGroups_", String(deg), "_",
                        String(nr), ".g.gz");
  strm := InputTextFile(path);
  r := EvalString(ReadAll(strm));
  CloseStream(strm);
  if r.degree <> deg or r.id <> nr then
    Error("record ", path, " is degree ", r.degree, " entry ", r.id);
  fi;
  if not IsBound(r.name) then
    r.name := "";
  fi;
  return [ r.id, r.size, r.SimpleSolvable, r.ONanScottType, r.suborbits,
           r.transitivity, r.name, r.SocleType, r.generators ];
end;

##  The vectors of GF(p)^d in the archive's numbering.
PRIMGRP_IntVectors := function(p, d)
  local vecs, i, r, v, j;
  vecs := [];
  for i in [1..p^d] do
    r := i - 1;
    v := [];
    for j in [1..d] do
      v[d-j+1] := r mod p;
      r := QuoInt(r, p);
    od;
    vecs[i] := v * Z(p)^0;
  od;
  return vecs;
end;

##  The permutation x -> x*m + t induces on the points, by arithmetic rather
##  than by searching the list of vectors.
PRIMGRP_AffinePerm := function(p, d, vecs, pw, m, t)
  local img, i, w, x, j;
  img := [];
  for i in [1..p^d] do
    w := vecs[i] * m + t;
    x := 1;
    for j in [1..d] do
      x := x + IntFFE(w[j]) * pw[j];
    od;
    img[i] := x;
  od;
  return PermList(img);
end;

##  The matrices of an affine group of degree p^d, or fail.  Point 1 is the
##  zero vector, so 1^g gives the translation part of a generator g, and the
##  images of the basis vectors, less that, give its rows.  Each generator must
##  be the affine map read off it, and the matrices must act irreducibly: then
##  they and one translation generate the group, as the reader assumes.
PRIMGRP_AffineMatrices := function(gens, deg, vecs, pw)
  local fac, p, d, one, basis, mats, g, t, m;
  fac := Factors(deg);
  p := fac[1];
  d := Length(fac);
  one := IdentityMat(d, GF(p));
  basis := List([1..d], j -> 1 + pw[j]);
  mats := [];
  for g in gens do
    t := vecs[1^g];
    m := List([1..d], j -> vecs[basis[j]^g] - t);
    if PRIMGRP_AffinePerm(p, d, vecs, pw, m, t) <> g then
      return fail;
    fi;
    if m <> one then
      Add(mats, m);
    fi;
  od;
  if not MTX.IsIrreducible(GModuleByMats(mats, d, GF(p))) then
    return fail;
  fi;
  return mats;
end;

PRIMGRP_Count := rec(alt := 0, sym := 0, prime := 0, matrices := 0,
                     generators := 0, named := 0);

##  One entry, as the line that stands in data/.
PRIMGRP_ImportEntry := function(deg, nr, vecs, pw)
  local l, c, ord, d, mats, p;
  l := PRIMGRP_ReadRecord(deg, nr);

  # The archive gives the order of the natural groups as a string such as
  # "Factorial(4097)/2".  Everything else the constructor gives must agree.
  if l[9] = "Alt" or l[9] = "Sym" then
    if l[9] = "Alt" then
      c := PGAlt(deg, nr);
      PRIMGRP_Count.alt := PRIMGRP_Count.alt + 1;
    else
      c := PGSym(deg, nr);
      PRIMGRP_Count.sym := PRIMGRP_Count.sym + 1;
    fi;
    ord := l[2];
    if IsString(ord) then
      ord := EvalString(ord);
    fi;
    if c[1] <> l[1] or c[2] <> ord or c{[3..6]} <> l{[3..6]}
       or c{[8,9]} <> l{[8,9]} then
      Error("degree ", deg, " entry ", nr, ": ", l[9],
            " does not agree with the archive");
    fi;
    if c[7] <> l[7] then
      if l[7] <> "" then
        Error("degree ", deg, " entry ", nr, ": named ", l[7], " not ", c[7]);
      fi;
      PRIMGRP_Count.named := PRIMGRP_Count.named + 1;
    fi;
    return Concatenation("[\"", l[9], "\"]");
  fi;

  if not IsInt(l[2]) then
    Error("degree ", deg, " entry ", nr, ": order ", l[2], " is not an integer");
  fi;

  if l[4] <> "1" then
    PRIMGRP_Count.generators := PRIMGRP_Count.generators + 1;
    return PRIMGRP_Compact(l);
  fi;

  if IsPrimeInt(deg) then
    d := l[2] / deg;
    c := PGPrime(d)(deg, nr);
    if c{[1..6]} <> l{[1..6]} or c[8] <> l[8] then
      Error("degree ", deg, " entry ", nr, ": [\"Prime\",", d,
            "] does not agree with the archive");
    fi;
    if c[7] <> l[7] then
      if l[7] <> "" then
        Error("degree ", deg, " entry ", nr, ": named ", l[7], " not ", c[7]);
      fi;
      PRIMGRP_Count.named := PRIMGRP_Count.named + 1;
    fi;
    PRIMGRP_Count.prime := PRIMGRP_Count.prime + 1;
    return Concatenation("[\"Prime\",", String(d), "]");
  fi;

  mats := PRIMGRP_AffineMatrices(l[9], deg, vecs, pw);
  if mats = fail or mats = [] then
    Error("degree ", deg, " entry ", nr, ": no affine matrices");
  fi;
  p := Factors(deg)[1];
  PRIMGRP_Count.matrices := PRIMGRP_Count.matrices + 1;
  return Concatenation(PRIMGRP_Compact(l{[1..8]}){[1..Length(PRIMGRP_Compact(l{[1..8]}))-1]},
           ",Z(", String(p), ")^0*",
           PRIMGRP_Compact(List(mats, m -> List(m, r -> List(r, IntFFE)))), "]");
end;

##  Write the degrees of one file; return the number of entries written.
PRIMGRP_ImportFile := function(path, degs)
  local out, n, deg, fac, vecs, pw, nr;
  out := OutputTextFile(path, false);
  SetPrintFormattingStatus(out, false);
  n := 0;
  for deg in degs do
    fac := Factors(deg);
    if IsPrimePowerInt(deg) and not IsPrimeInt(deg) then
      vecs := PRIMGRP_IntVectors(fac[1], Length(fac));
      pw := List([1..Length(fac)], j -> fac[1]^(Length(fac)-j));
    else
      vecs := fail;
      pw := fail;
    fi;
    PrintTo(out, "PRIMGRP[", deg, "]:=[\n");
    for nr in [1..NrPrimitiveGroups(deg)] do
      PrintTo(out, PRIMGRP_ImportEntry(deg, nr, vecs, pw), ",\n");
      n := n + 1;
    od;
    PrintTo(out, "];\n");
  od;
  CloseStream(out);
  return n;
end;

##  Append to PRIMINDX in <path> the file number of each degree in <indx>,
##  which the reader looks up before anything else: PrimGrpLoad refuses a
##  degree PRIMINDX does not bind, even one inside PRIMRANGE.  Re-running the
##  import over a file already extended is refused rather than doubling it.
PRIMGRP_ExtendIndex := function(path, indx)
  local s, i, j, k, old, rows, row, v;
  s := StringFile(path);
  if s = fail then
    Error("cannot read ", path);
  fi;
  i := PositionSublist(s, "BindGlobal(\"PRIMINDX\",");
  if i = fail then
    Error(path, " does not bind PRIMINDX");
  fi;
  j := Position(s, '[', i);
  k := PositionSublist(s, "]\n);", j);
  if k = fail then
    Error("the PRIMINDX list in ", path, " does not end as expected");
  fi;
  old := Filtered(SplitString(s{[j+1..k-1]}, ",\n"), x -> x <> "");
  if Length(old) <> 4095 then
    Error("PRIMINDX holds ", Length(old), " degrees, not the 4095 below the ",
          "import; extend it from the state before an earlier import");
  fi;

  # rows no wider than the ones already there
  rows := [];
  row := "";
  for v in indx do
    if Length(row) + Length(String(v)) >= 72 then
      Add(rows, row);
      row := "";
    fi;
    Append(row, String(v));
    Add(row, ',');
  od;
  Add(rows, row{[1..Length(row)-1]});
  FileString(path, Concatenation(s{[1..k-1]}, ",\n",
             JoinStringsWithSeparator(rows, "\n"), s{[k..Length(s)]}));
  Print("PRIMINDX: ", Length(old), " -> ", Length(old) + Length(indx),
        " degrees, in ", path, "\n");
end;

##  Split the degrees across files of about imp_perfile entries each -- a
##  degree is never split, since PRIMINDX maps it to one file -- write them,
##  and extend PRIMINDX with the file each degree went to.
PRIMGRP_ImportAll := function()
  local blocks, cur, size, deg, k, n, i, indx;
  blocks := [];
  cur := [];
  size := 0;
  for deg in [4096..8191] do
    k := NrPrimitiveGroups(deg);
    if size > 0 and size + k > imp_perfile then
      Add(blocks, cur);
      cur := [];
      size := 0;
    fi;
    Add(cur, deg);
    size := size + k;
  od;
  Add(blocks, cur);

  n := 0;
  indx := [];
  for i in [1..Length(blocks)] do
    n := n + PRIMGRP_ImportFile(
               Concatenation(imp_dir, "/gps", String(imp_first + i - 1), ".g"),
               blocks[i]);
    for deg in blocks[i] do
      indx[deg] := imp_first + i - 1;
    od;
    Print("FILE gps", imp_first + i - 1, " degrees ", blocks[i][1], " to ",
          blocks[i][Length(blocks[i])], ", ", n, " entries so far\n");
  od;
  PRIMGRP_ExtendIndex(imp_grp, indx{[4096..8191]});
  Print("IMPORTED ", n, " entries into ", Length(blocks), " files: ",
        PRIMGRP_Count, "\n");
end;

PRIMGRP_ImportAll();
QUIT;
