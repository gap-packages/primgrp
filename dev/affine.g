#############################################################################
##
##  dev/affine.g
##
##  Recover the compact matrix encoding for affine primitive groups that are
##  stored as explicit permutations.  Not read by read.g -- this is only used
##  to regenerate the data files.
##
##  An affine primitive group of degree n = p^d has O_p(G) as its socle: an
##  elementary abelian regular normal subgroup.  The library identifies point i
##  with Elements(GF(p)^d)[i], so the socle acts by translation and the point
##  stabiliser acts linearly.  Storing those matrices is enough; the reader
##  regenerates the translations from a single one (lib/primitiv.gi:121-134).
##

##  Two enumerations of F_p^d are in play, and they are not the same one:
##
##    "ffe" -- Elements(GF(p)^d), sorted by GAP's order on field elements,
##             which for a prime field is 0, Z(p)^0, Z(p)^1, Z(p)^2, ...
##             i.e. by discrete logarithm.  This is what lib/primitiv.gi:124
##             uses, so it is the convention of data/gps*.g.
##
##    "int" -- by the integer the base-p digits spell out: 0, 1, 2, ..., p-1.
##             This is the convention of the degree 4096-8191 data.
##
##  For p = 2 and p = 3 the two coincide, which is why degrees 4096 = 2^12 and
##  6561 = 3^8 -- between them 77% of the external data -- look the same either
##  way.  For p >= 5 they differ, and picking the wrong one silently yields a
##  conjugate group instead of the stored one.

BindGlobal("PRIMGRP_AffineContextCache", rec());

# PRIMGRP_AffineVectors is provided by lib/primitiv.gi, which documents the
# two enumerations; the conversion has to agree with the reader exactly, so it
# uses that one rather than a copy.

BindGlobal("PRIMGRP_AffineContext", function(deg, enum)
  local key, fac, p, d, vecs, pos, i, j, idx;
  key := Concatenation(String(deg), enum);
  if IsBound(PRIMGRP_AffineContextCache.(key)) then
    return PRIMGRP_AffineContextCache.(key);
  fi;
  fac := Factors(deg);
  p := fac[1];
  d := Length(fac);
  if not IsPrimeInt(p) or p^d <> deg then
    return fail;   # not a prime power, so not affine
  fi;
  vecs := PRIMGRP_AffineVectors(p, d, enum);
  # Position lookup by array, not by NewDictionary.  For a large prime field
  # the dictionary degrades to a linear scan, which turns the O(deg) labelling
  # check into O(deg^2): degree 7211 cost 386 s, degree 6561 cost 1 s.  Every
  # vector has a base-p index, so the position is just an array entry.
  pos := ListWithIdenticalEntries(deg, 0);
  for i in [1..deg] do
    idx := 0;
    for j in [1..d] do
      idx := idx * p + IntFFE(vecs[i][j]);
    od;
    pos[idx+1] := i;
  od;
  PRIMGRP_AffineContextCache.(key) := rec(p := p, d := d, enum := enum,
                                          vecs := vecs, pos := pos);
  return PRIMGRP_AffineContextCache.(key);
end);

# Which point is the vector <v>?
BindGlobal("PRIMGRP_VecPos", function(c, v)
  local idx, j;
  idx := 0;
  for j in [1..c.d] do
    idx := idx * c.p + IntFFE(v[j]);
  od;
  return c.pos[idx+1];
end);

# Reading off the affine structure.
#
# G is affine primitive of degree p^d, so under the right identification of
# the points with F_p^d every g in G is x -> x*A_g + w_g.  Both parts are
# computable straight from g: w_g is the image of zero, and A_g is read off
# from the images of the standard basis vectors.  Verifying that the formula
# holds at every point then proves G <= AGL(d,p) for this labelling.
#
# Nothing here needs a stabiliser chain of G, and that is the point: G has
# degree up to 8191 and order in the millions, and deterministic Schreier-Sims
# on it costs 44 s per group -- PCore, `t in G' and Stabilizer(G,1) each paid
# that toll, and between them they were the entire runtime of the conversion.
# The matrix group is d-dimensional and its order is immediate.
#
# Returns fail if <G> is not affine with this labelling.
BindGlobal("PRIMGRP_AffineMatrices", function(G, deg, enum)
  local c, basis, mats, g, w, m, ffm, i, j;
  c := PRIMGRP_AffineContext(deg, enum);
  if c = fail then return fail; fi;
  if not IsZero(c.vecs[1]) then return fail; fi;
  basis := List(IdentityMat(c.d, GF(c.p)), v -> PRIMGRP_VecPos(c, v));
  mats := [];
  for g in GeneratorsOfGroup(G) do
    w := c.vecs[1^g];                                  # image of zero
    ffm := List([1..c.d], j -> c.vecs[basis[j]^g] - w);
    for i in [1..deg] do
      if c.vecs[i^g] - w <> c.vecs[i] * ffm then
        return fail;                                   # not linear: wrong labelling
      fi;
    od;
    Add(mats, List(ffm, r -> List(r, IntFFE)));
  od;
  return rec(mats := mats);
end);

# Rebuild exactly as lib/primitiv.gi does, and insist on equality.
BindGlobal("PRIMGRP_AffineGroupFromMatrices", function(mats, deg, enum)
  local c, ffmats, perms, t;
  c := PRIMGRP_AffineContext(deg, enum);
  if Length(mats) = 0 then
    return Image(IsomorphismPermGroup(CyclicGroup(deg)));
  fi;
  ffmats := List(mats, m -> ImmutableMatrix(GF(c.p), Z(c.p)^0 * m));
  perms := List(ffmats, m -> Permutation(m, c.vecs, OnRight));
  t := First(c.vecs, v -> not IsZero(v));
  Add(perms, Permutation(t, c.vecs, \+));
  return Group(perms);
end);

# Which enumeration does <G> use, if either?  Whichever is returned has had
# every generator checked at every point.
BindGlobal("PRIMGRP_AffineEnumeration", function(G, deg)
  local order, enum;
  if deg > 4095 then
    order := [ "int", "ffe" ];
  else
    order := [ "ffe", "int" ];
  fi;
  for enum in order do
    if PRIMGRP_AffineMatrices(G, deg, enum) <> fail then
      return enum;
    fi;
  od;
  return fail;
end);

##  Verifying the conversion.
##
##  PRIMGRP_AffineMatrices has already checked, at every point, that each
##  generator of G is x -> x*A + w.  So G sits inside T : L, where T is the
##  full translation group and L = <A_g> is the image of G in GL(d,p), and
##  that is precisely the group the reader builds.  Two things remain:
##
##    (a) L acts irreducibly, so the reader's single stored translation
##        regenerates all of T;
##    (b) |L| * deg = |G|, which forces G to contain T and hence to be all of
##        T : L rather than a proper subgroup.
##
##  <storedSize> is the order from the data, not a recomputation -- a wrong
##  one makes the check fail and the entry keeps its permutations.
BindGlobal("PRIMGRP_CheckAffineRoundTrip", function(rec_, deg, enum, storedSize)
  local c, ffmats;
  c := PRIMGRP_AffineContext(deg, enum);
  if c = fail or Length(rec_.mats) = 0 then return false; fi;
  ffmats := List(rec_.mats, m -> ImmutableMatrix(GF(c.p), Z(c.p)^0 * m));
  if c.d > 1 and not MTX.IsIrreducible(GModuleByMats(ffmats, GF(c.p))) then
    return false;
  fi;
  return Size(Group(ffmats)) * deg = storedSize;
end);

# The expensive version, for spot-checking that the cheap one is honest.
BindGlobal("PRIMGRP_CheckAffineRoundTripSlow", function(G, deg, mats, enum)
  return PRIMGRP_AffineGroupFromMatrices(mats, deg, enum) = G;
end);
