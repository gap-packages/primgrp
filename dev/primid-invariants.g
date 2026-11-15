#############################################################################
##
##  This file is part of the PrimGrp package for GAP.
##
##  Copyright of PrimGrp belongs to its developers; see the manual,
##  PackageInfo.g and the git log for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
##  Compute data/primid.g: for the library groups that the fields stored in
##  data/ do not tell apart, data that PrimitiveIdentification checks against
##  random elements of the group to identify, so that it builds no library
##  group.
##
##  A class is a set of at least two groups of one degree that agree in
##  order, O'Nan-Scott type, suborbit lengths and transitivity, of a type other
##  than 4c, which is identified exactly.  For each class the data are
##
##      [ nrs, signs, types, masks ]
##
##    nrs     the numbers of its groups
##    signs   SignPermGroup of each
##    types   the cycle types, as PRIMGRP_CycleTypeHash, that occur in some
##            but not all of its groups
##    masks   for each group the types occurring in it, as an integer with
##            bit i-1 set for types[i]
##
##  An element of a conjugate of group j whose cycle type is among `types'
##  rules out every group lacking that type; the others are common to all.
##
##  The cycle types come from ConjugacyClasses, except for the classes
##  {PGL(2,q), M(q)} in degree q+1, q = r^2 odd, for which they follow from
##  the known classes of PGL(2,q) and PGL(2,r); ConjugacyClasses takes days
##  for q = 7921.
##
##  The work is split into one unit per class, so it can be spread over many
##  processes and resumed.  From the package directory, with GAP started as
##  gap -q -b -A --quitonbreak -l "ROOT;":
##
##    -c 'mode := "units";; out := "DIR";;' dev/primid-invariants.g
##        writes DIR/units.txt, one line "deg c" per class
##    -c 'mode := "unit";; out := "DIR";; deg := 7300;; c := 1;;' dev/primid-invariants.g
##        computes class c of degree deg into DIR/<deg>_<c>.g
##    -c 'mode := "degrees";; out := "DIR";; degs := [2..100];;' dev/primid-invariants.g
##        computes all classes of these degrees in one process, skipping
##        those already in DIR
##    -c 'mode := "assemble";; out := "DIR";;' dev/primid-invariants.g
##        checks that all units are there and writes data/primid.g
##
LoadPackage("primgrp");
SetInfoLevel(InfoWarning, 0);
SizeScreen([4096,]);

#############################################################################
##
##  Cycle types of PGL(2,q) and M(q) = PSL(2,q)<delta phi> on the q+1 points
##  of the projective line, q = r^2 odd, phi: x -> x^r, delta in PGL(2,q)
##  outside PSL(2,q).
##
##  An element of PGL(2,s) of order d > 1 is split, with 2 fixed points on
##  the line over GF(s), non-split, with none, or unipotent, with one; the
##  split orders are the divisors of s-1, the non-split ones those of s+1,
##  the unipotent one is the characteristic p, and in PSL(2,s) the
##  semisimple orders divide (s-1)/2 and (s+1)/2.
##
##  By Shintani descent the classes of PGL(2,q) phi correspond to those of
##  PGL(2,r), the coset PSL(2,q) delta phi to PGL(2,r) \ PSL(2,r): to
##  g phi belongs y in PGL(2,r) with (g phi)^2 conjugate to y.  An odd power
##  (g phi)^k fixes as many points of the line over GF(q) as y^k does of the
##  line over GF(r); the even power (g phi)^(2j) is conjugate to y^j, whose
##  semisimple elements are split over GF(q).
##

# the fixed points on the line over GF(s) of an element of order d and kind
PSL2Fixed := function(s, kind, d)
  if d = 1 then
    return s + 1;
  elif kind = "split" then
    return 2;
  elif kind = "nonsplit" then
    return 0;
  fi;
  return 1;
end;

# [ kind, order ] of the elements other than 1 of PGL(2,s), or of PSL(2,s)
PSL2Elements := function(s, p, inPSL)
  local m;
  m := 1;
  if inPSL then
    m := 2;
  fi;
  return Concatenation(
    List(Filtered(DivisorsInt((s-1)/m), d -> d > 1), d -> ["split", d]),
    List(Filtered(DivisorsInt((s+1)/m), d -> d > 1), d -> ["nonsplit", d]),
    [ ["unipotent", p] ]);
end;

# the cycle type, as CycleStructurePerm, of an element of order m whose k-th
# power fixes fix(k) points
CycleTypeFromFixed := function(m, fix)
  local cs, d, e;
  cs := [];
  for d in Filtered(DivisorsInt(m), d -> d > 1) do
    e := Sum(DivisorsInt(d), k -> MoebiusMu(d/k) * fix(k));    # points of period d
    if e mod d <> 0 then
      Error("points of period ", d, " do not form cycles");
    fi;
    if e > 0 then
      cs[d-1] := e / d;
    fi;
  od;
  return cs;
end;

# the cycle types of PGL(2,q) (m = false) or M(q) (m = true)
PSL2CycleTypes := function(q, m)
  local p, r, inner, outer;
  p := SmallestRootInt(q);
  r := RootInt(q);
  inner := List(PSL2Elements(q, p, m), x -> CycleTypeFromFixed(x[2],
             k -> PSL2Fixed(q, x[1], x[2] / Gcd(x[2], k))));
  outer := [];
  if m then
    outer := List(Difference(PSL2Elements(r, p, false), PSL2Elements(r, p, true)),
      y -> CycleTypeFromFixed(2 * y[2], function(k)
        if IsOddInt(k) then
          return PSL2Fixed(r, y[1], y[2] / Gcd(y[2], k));
        elif y[1] = "unipotent" then
          return PSL2Fixed(q, y[1], y[2] / Gcd(y[2], k/2));
        fi;
        return PSL2Fixed(q, "split", y[2] / Gcd(y[2], k/2));
      end));
  fi;
  return Set(Concatenation([ [] ], inner, outer));
end;

# the data of class nrs of degree deg
ClassData := function(deg, nrs)
  local q, e, pgl;
  q := deg - 1;
  e := List(nrs, i -> PRIMGrp(deg, i));
  pgl := Concatenation("PGL(2,", String(q), ")");
  if Length(nrs) = 2 and IsOddInt(q) and IsPrimePowerInt(q)
     and RootInt(q)^2 = q and pgl in List(e, x -> x[7]) and e[1][6] = 3 then
    return PRIMGRP_IdClassDataFromTypes(nrs,
             List(nrs, i -> SignPermGroup(PrimitiveGroup(deg, i))),
             List(e, x -> PSL2CycleTypes(q, x[7] <> pgl)));
  fi;
  return PRIMGRP_IdClassData(deg, nrs);
end;

if mode = "units" then
  PrintTo(Concatenation(out, "/units.txt"));
  for deg in [2..8191] do
    for c in [1..Length(PRIMGRP_IdClasses(deg))] do
      AppendTo(Concatenation(out, "/units.txt"), deg, " ", c, "\n");
    od;
  od;

elif mode = "unit" then
  nrs := PRIMGRP_IdClasses(deg)[c];
  t := Runtime();
  data := ClassData(deg, nrs);
  # written to a temporary file first, so a killed run leaves no unit behind
  file := Concatenation(out, "/", String(deg), "_", String(c), ".g");
  PrintTo(Concatenation(file, ".tmp"), "PRIMGRP_IDUNIT := ", [deg, c, data], ";\n",
          "PRIMGRP_IDTIME := ", Runtime() - t, ";\n");
  Exec(Concatenation("mv ", file, ".tmp ", file));

elif mode = "degrees" then
  for deg in degs do
    for c in [1..Length(PRIMGRP_IdClasses(deg))] do
      file := Concatenation(out, "/", String(deg), "_", String(c), ".g");
      if not IsExistingFile(file) then
        PrintTo(Concatenation(file, ".tmp"), "PRIMGRP_IDUNIT := ",
                [deg, c, ClassData(deg, PRIMGRP_IdClasses(deg)[c])], ";\n");
        Exec(Concatenation("mv ", file, ".tmp ", file));
      fi;
    od;
  od;

elif mode = "assemble" then
  units := List(SplitString(StringFile(Concatenation(out, "/units.txt")), "\n"),
                l -> List(SplitString(l, " "), Int));
  if IsBound(maxdeg) then
    units := Filtered(units, u -> u[1] <= maxdeg);    # for testing
  fi;
  res := [];
  for u in units do
    Read(Concatenation(out, "/", String(u[1]), "_", String(u[2]), ".g"));
    if PRIMGRP_IDUNIT{[1,2]} <> u then
      Error("unit file ", u, " holds ", PRIMGRP_IDUNIT{[1,2]});
    fi;
    if not IsBound(res[u[1]]) then
      res[u[1]] := [];
    fi;
    res[u[1]][u[2]] := PRIMGRP_IDUNIT[3];
  od;
  file := Filename(DirectoriesPackageLibrary("primgrp", "data")[1], "primid.g");
  PrintTo(file, "# written by dev/primid-invariants.g\n");
  for deg in [1..Length(res)] do
    if IsBound(res[deg]) then
      AppendTo(file, "PRIMGRP_IDDATA.data[", deg, "]:=",
               ReplacedString(String(res[deg]), " ", ""), ";\n");
    fi;
  od;
fi;
QUIT;
