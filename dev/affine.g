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
  local key, fac, p, d, vecs, pos, i;
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
  pos := NewDictionary(vecs[1], true, vecs);
  for i in [1..deg] do
    AddDictionary(pos, vecs[i], i);
  od;
  PRIMGRP_AffineContextCache.(key) := rec(p := p, d := d, enum := enum,
                                          vecs := vecs, pos := pos);
  return PRIMGRP_AffineContextCache.(key);
end);

# Does the socle act as translation by v[i] under this enumeration?  This is
# what makes the conversion lossless rather than merely correct up to conjugacy.
BindGlobal("PRIMGRP_HasLabelling", function(G, deg, enum)
  local c, T, t, shift, i;
  c := PRIMGRP_AffineContext(deg, enum);
  if c = fail then return false; fi;
  T := PCore(G, c.p);
  if Size(T) <> deg then return false; fi;
  if not IsZero(c.vecs[1]) then return false; fi;
  for t in GeneratorsOfGroup(T) do
    shift := c.vecs[1^t];
    for i in [1..deg] do
      if i^t <> LookupDictionary(c.pos, c.vecs[i] + shift) then
        return false;
      fi;
    od;
  od;
  return true;
end);

# Which enumeration does <G> use, if either?
BindGlobal("PRIMGRP_AffineEnumeration", function(G, deg)
  local enum;
  for enum in [ "ffe", "int" ] do
    if PRIMGRP_HasLabelling(G, deg, enum) then
      return enum;
    fi;
  od;
  return fail;
end);

# The matrices of a small generating set of the point stabiliser, as integer
# matrices (the reader scales them by Z(p)^0).  Returns fail if <G> is not
# affine with this labelling.
#
# Returns the generating set alongside the matrices.  SmallGeneratingSet is
# randomised, so a second call gives a different set and the matrices would no
# longer correspond to it -- the two have to travel together.
BindGlobal("PRIMGRP_AffineMatrices", function(G, deg, enum)
  local c, S, basis, mats, s, m, j;
  c := PRIMGRP_AffineContext(deg, enum);
  if c = fail then return fail; fi;
  if not PRIMGRP_HasLabelling(G, deg, enum) then return fail; fi;
  # The stabiliser can be trivial (G is then the regular group C_p), so keep
  # the generating list rather than building a group from it.
  S := SmallGeneratingSet(Stabilizer(G, 1));
  # point number of the j-th standard basis vector
  basis := List(IdentityMat(c.d, GF(c.p)), v -> LookupDictionary(c.pos, v));
  mats := [];
  for s in S do
    m := List([1..c.d], j -> List(c.vecs[basis[j]^s], IntFFE));
    Add(mats, m);
  od;
  return rec(gens := S, mats := mats);
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

##  Verifying the conversion.
##
##  Rebuilding the degree 8191 permutation group and testing `= G' needs two
##  stabiliser chains and costs seconds per group; over 19571 groups that is
##  days.  The same statement can be settled in the d-dimensional matrix group,
##  which is tiny.  If
##
##    (a) the socle acts as translation by v[i]        [PRIMGRP_HasLabelling]
##    (b) each stabiliser generator s satisfies v[i^s] = v[i]*M_s for every i
##    (c) <M_s> acts irreducibly on F_p^d
##    (d) |<M_s>| * deg = |G|
##
##  then G = T : G_1 = <translations, matrices> is exactly what the reader
##  builds: (b) says the matrices are the linear action, (c) says the single
##  stored translation generates all of T under it, and (a) and (d) say nothing
##  is missing.  <storedSize> is the size from the data, not a recomputation.
BindGlobal("PRIMGRP_CheckAffineRoundTrip", function(rec_, deg, enum, storedSize)
  local c, S, mats, ffmats, i, j, s, m;
  c := PRIMGRP_AffineContext(deg, enum);
  if c = fail then return false; fi;
  S := rec_.gens;
  mats := rec_.mats;
  ffmats := List(mats, m -> ImmutableMatrix(GF(c.p), Z(c.p)^0 * m));
  for j in [1..Length(S)] do                                   # (b)
    s := S[j]; m := ffmats[j];
    for i in [1..deg] do
      if i^s <> LookupDictionary(c.pos, c.vecs[i] * m) then
        return false;
      fi;
    od;
  od;
  if Length(mats) = 0 then                                     # G is C_p
    return storedSize = deg;
  fi;
  if not MTX.IsIrreducible(GModuleByMats(ffmats, GF(c.p))) then # (c)
    return false;
  fi;
  return Size(Group(ffmats)) * deg = storedSize;               # (d)
end);

# The expensive version, for spot-checking that the cheap one is honest.
BindGlobal("PRIMGRP_CheckAffineRoundTripSlow", function(G, deg, mats, enum)
  return PRIMGRP_AffineGroupFromMatrices(mats, deg, enum) = G;
end);
