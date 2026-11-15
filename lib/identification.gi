#############################################################################
##
##  This file is part of the PrimGrp package for GAP.
##
##  Copyright of PrimGrp belongs to its developers; see the manual,
##  PackageInfo.g and the git log for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
##  PrimitiveIdentification and SimsNo.
##

# local cache for `PrimitiveIdentification':
BindGlobal("PGIdCache", rec( deg := 0, data := [] ));

#############################################################################
##
#F  PRIMGRP_HasJordanElement( <G>, <dom> )
##
##  For a transitive group <G> on <dom>, true if a random element shows that
##  <G> contains Alt(n), n = |dom|; false otherwise.
##
##  An element with a cycle of prime length p, n/2 < p < n-2, has a power that
##  is a p-cycle, and makes <G> primitive; so <G> contains Alt(n) by Jordan's
##  theorem.  In Alt(n) and Sym(n), the cycle through a given point has such a
##  length with probability about 1/(2 log n).  For n = 7526, 50 log2(n) tries
##  miss with probability 10^-15; the first test in IsNaturalSymmetricGroup
##  tries 10 log2(n) times, misses once in a thousand, and then spends 78 s
##  in IsPrimitive.
##
BindGlobal("PRIMGRP_HasJordanElement", function(G, dom)
  local n, pnt, i, l;
  n := Length(dom);
  pnt := dom[1];
  for i in [1..50*LogInt(n, 2)] do
    l := CycleLength(PseudoRandom(G), pnt);
    if l > n/2 and l < n-2 and IsPrimeInt(l) then
      return true;
    fi;
  od;
  return false;
end);

InstallMethod(PrimitiveIdentification,"generic",true,[IsPermGroup],0,
function(grp)
local dom,deg,PD,s,cand,a,p,b,f,q,gl,hom;
  dom:=MovedPoints(grp);
  if Length(dom) < 2 then
    Error("PrimitiveIdentification: <G> must move at least two points");
  fi;
  deg:=Length(dom);

  # treat alternating and symmetric group: these are uniquely determined by
  # their size, and primitive.  Test them before IsPrimitive, which costs far
  # more for them: 36 s against 1 ms for Sym(3894).
  # IsNaturalAlternatingGroup and IsNaturalSymmetricGroup are not used here:
  # for other groups they go on to IsPrimitive and Size, which the code below
  # avoids for solvable groups, and until gap-system/gap#6850 their first
  # random search misses Alt(n) or Sym(n) once in about 1000 calls and then
  # runs IsPrimitive too.
  if not HasSize(grp) and IsTransitive(grp, dom)
     and PRIMGRP_HasJordanElement(grp, dom) then
    if ForAll(GeneratorsOfGroup(grp), x -> SignPerm(x) = 1) then
      s := Factorial(deg)/2;
      SetIsNaturalAlternatingGroup(grp, true);
    else
      s := Factorial(deg);
      SetIsNaturalSymmetricGroup(grp, true);
    fi;
    SetSize(grp, s);
    return First([1..NrPrimitiveGroups(deg)], i -> PRIMGrp(deg, i)[2] = s);
  fi;

  # For other groups the test below ends in Size, whose Schreier-Sims
  # verification takes 143 s for D(2*8191); a solvable group gets it from a
  # pcgs in 0.1 s instead, and IsSolvableGroup costs at most 55 ms otherwise.
  if not HasSize(grp) and not IsSolvableGroup(grp) then
    if IsNaturalAlternatingGroup(grp) then
      SetSize(grp, Factorial(deg)/2);
    elif IsNaturalSymmetricGroup(grp) then
      SetSize(grp, Factorial(deg));
    fi;
  fi;

  if not (IsTransitive(grp,dom) and IsPrimitive(grp,dom)) then
    Error("Group must operate primitively");
  fi;

  # get information of all primitive groups of the given degree
  PD := List([1 .. NrPrimitiveGroups(deg)], t -> PRIMGrp(deg, t));

  #
  # filter by size
  #
  s := Size(grp);
  PD := Filtered(PD, primdata -> primdata[2] = s);
  if Length(PD)=1 then return PD[1][1]; fi;

  #
  # filter by ONanScottType if it is not unique
  #
  if Length(Set(PD, primdata -> primdata[4])) > 1 then
    a := ONanScottType(grp);
    PD := Filtered(PD, primdata -> primdata[4] = a);
  fi;
  if Length(PD)=1 then return PD[1][1]; fi;

  #
  # filter by suborbit lengths
  #
  if Length(Set(PD, primdata -> primdata[5])) > 1 then
    s := Stabilizer(grp, dom[1]);
    a := Set(Collected(OrbitLengths(s, dom{[2..Length(dom)]})));
    PD := Filtered(PD, primdata -> Set(primdata[5]) = a);
  fi;
  if Length(PD)=1 then return PD[1][1]; fi;

  #
  # filter by transitivity if it is not unique
  #
  if Length(Set(PD, primdata -> primdata[6])) > 1 then
    a := Transitivity(grp, dom);
    PD := Filtered(PD, primdata -> primdata[6] = a);
  fi;
  if Length(PD)=1 then return PD[1][1]; fi;

  # now we need to create the groups
  cand := List(PD, primdata -> PrimitiveGroup(deg, primdata[1]));

  # in product action case, some tests on the socle quotient.
  if ONanScottType(grp) = "4c" then
    # first we just identify its isomorphism type
    f := g -> IdGroup(FactorGroup(g, Socle(g)));
    a := f(grp);
    cand := Filtered(cand, g -> f(g) = a);
  fi;
  if Length(cand) = 1 then return PrimitiveIdentification(cand[1]); fi;

  #
  # AbelianInvariants
  #
  if Length(cand) > 1 then
    a := AbelianInvariants(grp);
    cand := Filtered(cand, g -> AbelianInvariants(g) = a);
  fi;

  #
  # sylow orbits
  #
  if Length(cand) > 1 then
    # sylow orbits
    gl := ShallowCopy(PrimeDivisors(Size(grp)));
    f := function(g)
      local s;
      s := SylowSubgroup(g, p);
      return Collected(OrbitLengths(s, MovedPoints(g)));
    end;
    while Length(cand) > 1 and Length(gl) > 0 do
      p := Remove(gl);   # largest remaining prime divisor
      a := f(grp);
      cand := Filtered(cand, g -> f(g) = a);
    od;
  fi;

  #
  # Some further tests for the sylow subgroups
  #
  if Length(cand) > 1 then
    for q in PrimeDivisors(Size(grp)/Size(Socle(grp))) do
      if q=1 then
        q:=2;
      fi;

      # central series
      f := function(g)
        local ag;
        ag:=Image(IsomorphismPcGroup(SylowSubgroup(g,q)));
        return List(LowerCentralSeries(ag), Size);
      end;
      a := f(grp);
      cand := Filtered(cand, g -> f(g) = a);

      if Length(cand) > 1 then
        # Frattini subgroup
        f := function(g)
          local ag;
          ag:=Image(IsomorphismPcGroup(SylowSubgroup(g,q)));
          return Size(FrattiniSubgroup(ag));
        end;
        a := f(grp);
        cand := Filtered(cand, g -> f(g) = a);
      fi;

      if Length(cand) > 1 and Size(SylowSubgroup(grp,q)) < 512 then
        # Isomorphism type of Sylow subgroup
        f := function(g)
          local ag;
          ag:=Image(IsomorphismPcGroup(SylowSubgroup(g,q)));
          return IdGroup(ag);
        end;
        a := f(grp);
        cand := Filtered(cand, g -> f(g) = a);
      fi;
    od;
  fi;

  # back for a closer look at the product action groups.
  if Length(cand) > 1 and ONanScottType(grp) = "4c" then
    # construct the action of the socle quotient as a (necessarily transitive)
    # action on the socle factors.
    f := function(g)
      local s, cs, n, beta, alpha;
      s := Socle(g);
      cs := CompositionSeries(s);
      cs := cs[Length(cs)-1];
      n := Normalizer(g, cs);
      beta := FactorCosetAction(g, n);
      alpha := FactorCosetAction(n, ClosureGroup(Centralizer(n, cs), s));
      return TransitiveIdentification(Group(KuKGenerators(g, beta, alpha)));
    end;
    a := f(grp);
    cand := Filtered(cand, g -> f(g) = a);
  fi;

  #
  # conjugacy classes
  #
  if Length(cand) > 1 then

    # use caching
    if deg <> PGIdCache.deg then
      PGIdCache.deg := deg;
      PGIdCache.data := [];
    fi;

    # compute the fingerprint, with cache lookup
    f := function(g)
      local i, dat;
      if HasPrimitiveIdentification(g) then
        # is it in the cache?
        i := PrimitiveIdentification(g);
        if IsBound(PGIdCache.data[i]) then
          return PGIdCache.data[i];
        fi;
      fi;
      dat := ConjugacyClasses(g : onlysizes);
      dat := List(dat, i -> [CycleStructurePerm(Representative(i)), Size(i)]);
      dat := Collected(dat);
      if IsBound(i) then
        PGIdCache.data[i] := dat;
      fi;
      return dat;
    end;

    a := f(grp);
    cand := Filtered(cand, g -> f(g) = a);
  fi;

  #
  # affine case
  #
  if Length(cand) > 1 and ONanScottType(grp) = "1" then
    Assert(0, ForAll(cand, g -> ONanScottType(g) = "1"));

    # degree must be a prime power
    gl := Factors(deg);
    gl := GL(Length(gl), gl[1]);
    hom := IsomorphismPermGroup(gl);  # nasty

    # map group by its linear action into GL, and then further via hom
    f := g -> Image(hom, Subgroup(gl, LinearActionLayer(g, Pcgs(Socle(g)))));
    a := f(grp);
    b := Image(hom, gl);
    cand := Filtered(cand, g -> RepresentativeAction(b, f(g), a) <> fail);
  fi;

  if Length(cand) > 1 then
    Error("Uh-Oh, this should never happen ",cand);
  fi;
  return PrimitiveIdentification(cand[1]);
end);

InstallMethod(SimsNo,"via `PrimitiveIdentification'",true,[IsPermGroup],0,
function(grp)
local dom;
  dom:=MovedPoints(grp);
  if NrMovedPoints(grp) > 50 then
    Error("SimsNo is defined only for primitive groups of degree <= 50");
  fi;
  if not (IsTransitive(grp,dom) and IsPrimitive(grp,dom)) then
    Error("Group must operate primitively");
  fi;
  return SimsNo(PrimitiveGroup(Length(dom),PrimitiveIdentification(grp)));
end);
