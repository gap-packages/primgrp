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

#############################################################################
##
##  Groups of type 4c, that is, in product action on the m^k tuples over a set
##  of size m, with socle T^k for a non-abelian simple T that acts on each
##  coordinate.  Such a group lies in W = Sym(m) wr Sym(k), and two of them are
##  conjugate in Sym(m^k) if and only if they are conjugate in W: a conjugating
##  permutation maps socle factors to socle factors, hence the fibres of the
##  coordinates to each other, and that makes it an element of W.
##
##  W also acts on the m*k pairs (i,c) of a coordinate and a value, and there
##  a group of order 10^100 is no harder than its degree suggests.  All of the
##  following happens on the pairs; the library stores these groups by their
##  elements (p_1,...,p_k; sigma) of W in the first place.
##

#############################################################################
##
#F  PRIMGRP_WreathToPerm( <m>, <k>, <e> )
#F  PRIMGRP_PermToWreath( <m>, <k>, <w> )
##
##  Turn the element <e> = [ p_1, ..., p_k, sigma ] of Sym(m) wr Sym(k) into
##  the permutation of the pairs it induces, the pair (i,c) being the point
##  (i-1)*m+c; and back.  As in PGProductAction4c, <e> sends (i,c) to
##  (i^sigma, c^p_i).
##
BindGlobal("PRIMGRP_WreathToPerm", function(m, k, e)
  local img, i, c;
  img := [];
  for i in [1..k] do
    for c in [1..m] do
      img[(i-1)*m+c] := (i^e[k+1]-1)*m + c^e[i];
    od;
  od;
  return PermList(img);
end);

BindGlobal("PRIMGRP_PermToWreath", function(m, k, w)
  local e, sigma, i, j;
  e := [];
  sigma := [];
  for i in [1..k] do
    j := QuoInt(((i-1)*m+1)^w - 1, m) + 1;
    sigma[i] := j;
    e[i] := PermList(List([1..m], c -> ((i-1)*m+c)^w - (j-1)*m));
  od;
  e[k+1] := PermList(sigma);
  return e;
end);

#############################################################################
##
#F  PRIMGRP_CoordinateFibres( <G>, <m>, <k>, <e> )
##
##  For <G> of type 4c on m^k points, return for each coordinate the partition
##  of the points by their value there, or fail.  <e> must be a multiple of the
##  index of the socle.
##
##  Only orbits of groups given by a few generators are computed, never a
##  stabiliser chain, which for a socle like Alt(69)^2 is the difference
##  between milliseconds and minutes.  The price is that the answer may be
##  wrong: the generators taken need not generate what they are meant to.
##  PRIMGRP_WreathElements checks it.
##
BindGlobal("PRIMGRP_CoordinateFibres", function(G, m, k, e)
  local deg, x, S, y, o, U, lines, comps, parts, tries, g, V, q, fib, i,
        others;
  deg := m^k;

  # G/Soc(G) has order dividing e, so an e-th power lies in the socle; and
  # the socle is the normal closure of any of its non-trivial elements
  repeat
    x := PseudoRandom(G)^e;
  until not IsOne(x);
  S := Group(List([1..12], i -> x^PseudoRandom(G)));

  # an element of prime order has a fair chance to lie in a single factor
  # T_i.  Its conjugates under the socle then span a subgroup of T_i, whose
  # orbits are the lines in direction i.
  y := PseudoRandom(S);
  o := Order(y);
  if o = 1 then
    return fail;
  fi;
  y := y^(o/Random(PrimeDivisors(o)));
  U := Group(List([1..10], i -> y^PseudoRandom(S)));
  lines := Orbits(U, [1..deg]);
  if Length(lines) <> m^(k-1) or ForAny(lines, l -> Length(l) <> m) then
    return fail;
  fi;

  # G permutes the factors transitively; conjugates of U are told apart by
  # their lines
  comps := [U];
  parts := [Set(lines, Set)];
  tries := 0;
  while Length(comps) < k and tries < 100 do
    tries := tries + 1;
    g := PseudoRandom(G);
    V := Group(List(GeneratorsOfGroup(U), u -> u^g));
    q := Set(Orbits(V, [1..deg]), Set);
    if not q in parts then
      Add(comps, V);
      Add(parts, q);
    fi;
  od;
  if Length(comps) < k then
    return fail;
  fi;

  # the fibres of coordinate i are the orbits of all factors but the i-th
  fib := [];
  for i in [1..k] do
    others := Group(Concatenation(List(Filtered([1..k], j -> j <> i),
                                       j -> GeneratorsOfGroup(comps[j]))));
    fib[i] := Set(Orbits(others, [1..deg]), Set);
    if Length(fib[i]) <> m or ForAny(fib[i], f -> Length(f) <> m^(k-1)) then
      return fail;
    fi;
  od;
  return fib;
end);

#############################################################################
##
#F  PRIMGRP_WreathElements( <G>, <m>, <k>, <e> )
##
##  For <G> of type 4c on m^k points, return its generators as elements
##  [ p_1, ..., p_k, sigma ] of Sym(m) wr Sym(k), or fail if the random choices
##  in PRIMGRP_CoordinateFibres went wrong; then ask again.  <e> is as there.
##
##  What is returned has been verified: PGProductAction4c turns it back into
##  the generators of <G>, relabelled.
##
BindGlobal("PRIMGRP_WreathElements", function(G, m, k, e)
  local deg, fib, coord, i, a, x, pw, rel, els, g, sigma, ps, img, j;
  deg := m^k;
  fib := PRIMGRP_CoordinateFibres(G, m, k, e);
  if fib = fail then
    return fail;
  fi;

  # the coordinates of each point, and the relabelling that numbers the
  # points as PGProductAction4c does
  coord := List([1..deg], x -> []);
  for i in [1..k] do
    for a in [1..m] do
      for x in fib[i][a] do
        coord[x][i] := a;
      od;
    od;
  od;
  pw := List([1..k], i -> m^(k-i));
  rel := List(coord, c -> 1 + Sum([1..k], i -> (c[i]-1)*pw[i]));
  if Length(Set(rel)) <> deg then
    return fail;
  fi;
  rel := PermList(rel);

  # a generator sends the fibres of coordinate i to those of coordinate
  # j = i^sigma, and what it does to them is p_i
  els := [];
  for g in GeneratorsOfGroup(G) do
    sigma := [];
    ps := [];
    for i in [1..k] do
      img := OnTuples(fib[i][1], g);
      j := First([1..k], j -> ForAll(img, x -> coord[x][j] = coord[img[1]][j]));
      if j = fail then
        return fail;
      fi;
      sigma[i] := j;
      ps[i] := List([1..m], a -> coord[fib[i][a][1]^g][j]);
    od;
    if Length(Set(sigma)) <> k or ForAny(ps, p -> Length(Set(p)) <> m) then
      return fail;
    fi;
    Add(els, Concatenation(List(ps, PermList), [PermList(sigma)]));
  od;

  if GeneratorsOfGroup(PGProductAction4c(m, k, els))
     <> List(GeneratorsOfGroup(G), x -> x^rel) then
    return fail;
  fi;
  return els;
end);

#############################################################################
##
#F  PRIMGRP_AlignSocleFactors( <m>, <k>, <G> )
##
##  <G> is a group of type 4c acting on the m*k pairs.  Its socle is
##  T_1 x ... x T_k with T_i acting on the values of coordinate i.  The T_i are
##  conjugate in Sym(m), but not equal, since how the values of each
##  coordinate are numbered is arbitrary.
##
##  Return a record with a conjugate `G' of <G> under the base group Sym(m)^k
##  of the wreath product in which all T_i are equal to T_1, and with that
##  group as `T'.
##
BindGlobal("PRIMGRP_AlignSocleFactors", function(m, k, G)
  local H, T, cs, blocks, i, g;

  # what the stabiliser of the first coordinate does to its values
  H := Action(Stabilizer(G, [1..m], OnSets), [1..m]);
  T := Socle(H);

  # if g sends coordinate 1 to coordinate i, acting as p on its values, then
  # T_i = T_1^p
  cs := [()];
  blocks := ActionHomomorphism(G, List([1..k], i -> [(i-1)*m+1..i*m]), OnSets);
  for i in [2..k] do
    g := PreImagesRepresentative(blocks,
           RepresentativeAction(Image(blocks), 1, i));
    cs[i] := PRIMGRP_PermToWreath(m, k, g)[1]^-1;
  od;

  g := PRIMGRP_WreathToPerm(m, k, Concatenation(cs, [()]));
  return rec( G := Group(List(GeneratorsOfGroup(G), x -> x^g)), T := T );
end);

#############################################################################
##
#F  PRIMGRP_ConjugatingPermutation( <T1>, <T2> )
##
##  For transitive groups <T1> and <T2> on [ 1 .. m ], return a permutation c
##  with T1^c = T2, or fail if there is none.
##
##  Conjugation by c is an isomorphism that maps point stabilisers to point
##  stabilisers, and any isomorphism that does so comes from a permutation.
##  So take one isomorphism and vary it by outer automorphisms; inner ones map
##  a point stabiliser to another.  RepresentativeAction in Sym(m) finds the
##  same, but for two copies of PSL(3,5) on 31 points not within minutes.
##
BindGlobal("PRIMGRP_ConjugatingPermutation", function(T1, T2)
  local iso, stab, gens1, gens2, aut, alpha, U, b, c, queue, j, i, x;
  iso := IsomorphismGroups(T1, T2);
  if iso = fail then
    return fail;
  fi;
  stab := Image(iso, Stabilizer(T1, 1));
  gens1 := GeneratorsOfGroup(T1);
  aut := AutomorphismGroup(T2);
  for alpha in RightTransversal(aut, InnerAutomorphismsAutomorphismGroup(aut)) do
    # a subgroup of the order of a point stabiliser that fixes a point is one
    U := Image(alpha, stab);
    b := First(MovedPoints(T2),
               b -> ForAll(GeneratorsOfGroup(U), u -> b^u = b));
    if b <> fail then
      # c maps 1^g to b^(image of g)
      gens2 := List(gens1, g -> Image(alpha, Image(iso, g)));
      c := [];
      c[1] := b;
      queue := [1];
      for j in queue do
        for i in [1..Length(gens1)] do
          x := j^gens1[i];
          if not IsBound(c[x]) then
            c[x] := c[j]^gens2[i];
            Add(queue, x);
          fi;
        od;
      od;
      return PermList(c);
    fi;
  od;
  return fail;
end);

#############################################################################
##
#F  PRIMGRP_IdentifyProductAction( <G>, <PD> )
##
##  <G> is a primitive group of type 4c and <PD> a list of library entries of
##  that type and the degree and order of <G>.  Return those to which <G> is
##  conjugate in the symmetric group.
##
##  Once the socle factors of two such groups are aligned, both contain the
##  same T^k and lie in its normaliser N wr Sym(k) in W, where N is the
##  normaliser of T in Sym(m).  So they are conjugate in W if and only if their
##  images are conjugate in the quotient (N/T) wr Sym(k).  N/T embeds into
##  Out(T), so that group is tiny.
##
BindGlobal("PRIMGRP_IdentifyProductAction", function(G, PD)
  local deg, k, m, natural, align, libs, e, tries, els, inp, N, hom, nu, O,
        mq, Q, quotient, a;
  deg := NrMovedPoints(G);
  k := PD[1][8][3];
  m := RootInt(deg, k);

  # Alt(m) on m points is normal in Sym(m), so there is nothing to align, and
  # nothing to compute to know T.  With Alt(69) that saves seconds per group.
  natural := ForAll(PD, primdata -> primdata[8][1] = "A"
                                    and primdata[8][2] = m);
  align := function(els)
    local H;
    H := Group(List(els, x -> PRIMGRP_WreathToPerm(m, k, x)));
    if natural then
      return rec( G := H, T := AlternatingGroup(m) );
    fi;
    return PRIMGRP_AlignSocleFactors(m, k, H);
  end;
  libs := List(PD, primdata -> align(primdata[9][2]));

  # a multiple of the index of the socle.  The candidates need not agree on
  # their socle.
  e := Lcm(List(libs, r -> Size(G) / Size(r.T)^k));

  tries := 0;
  repeat
    tries := tries + 1;
    if tries > 100 then
      Error("cannot write <G> as a group in product action");
    fi;
    els := PRIMGRP_WreathElements(G, m, k, e);
  until els <> fail;
  inp := align(els);

  # the map from N wr Sym(k) onto (N/T) wr Sym(k), both acting on pairs
  if natural then
    O := Group((1,2));
    nu := function(p)
      if SignPerm(p) = 1 then
        return ();
      fi;
      return (1,2);
    end;
  else
    N := Normalizer(SymmetricGroup(m), inp.T);
    hom := NaturalHomomorphismByNormalSubgroup(N, inp.T);
    hom := hom * IsomorphismPermGroup(Image(hom));
    O := Image(hom);
    nu := p -> Image(hom, p);
  fi;
  mq := Maximum(1, LargestMovedPoint(O));
  Q := Group(Concatenation(
         List(GeneratorsOfGroup(O), o -> PRIMGRP_WreathToPerm(mq, k,
              Concatenation([o], List([2..k], i -> ()), [()]))),
         List(GeneratorsOfGroup(SymmetricGroup(k)), s ->
              PRIMGRP_WreathToPerm(mq, k,
              Concatenation(List([1..k], i -> ()), [s])))), ());
  quotient := function(gens)
    return Group(List(gens, function(w)
      w := PRIMGRP_PermToWreath(m, k, w);
      return PRIMGRP_WreathToPerm(mq, k,
               Concatenation(List([1..k], i -> nu(w[i])), [w[k+1]]));
    end), ());
  end;

  # the kernel of G -> Q is the intersection of G with T^k, so this checks
  # that the socle of <G> is what the above took it to be
  a := quotient(GeneratorsOfGroup(inp.G));
  if Size(a) * Size(inp.T)^k <> Size(G) then
    Error("the socle of <G> is not T^k");
  fi;

  return PD{Filtered([1..Length(PD)], function(i)
    local c;
    if libs[i].T = inp.T then
      c := ();
    else
      c := PRIMGRP_ConjugatingPermutation(libs[i].T, inp.T);
      if c = fail then
        return false;
      fi;
    fi;
    c := PRIMGRP_WreathToPerm(m, k, Concatenation(List([1..k], j -> c), [()]));
    return RepresentativeAction(Q, a,
             quotient(List(GeneratorsOfGroup(libs[i].G), x -> x^c))) <> fail;
  end)};
end);

#############################################################################
##
#F  PRIMGRP_CycleTypeHash( <cs> )
##
##  An integer standing for the cycle type <cs>, as CycleStructurePerm returns
##  it: cs[i] is the number of cycles of length i+1.  The data in
##  data/primid.g store cycle types this way; PRIMGRP_IdClassDataFromTypes
##  checks that no two cycle types in one class get the same value.
##
BindGlobal("PRIMGRP_CycleTypeHash", function(cs)
  local h, i;
  h := 1;
  for i in [1..Length(cs)] do
    if IsBound(cs[i]) then
      h := (h * 1000003 + i) mod 2147483647;
      h := (h * 1000003 + cs[i]) mod 2147483647;
    fi;
  od;
  return h;
end);

#############################################################################
##
#F  PRIMGRP_IdClasses( <deg> )
#F  PRIMGRP_IdClassData( <deg>, <nrs> )
#F  PRIMGRP_IdClassDataFromTypes( <nrs>, <signs>, <cts> )
##
##  The classes of degree <deg> that data/primid.g covers, as lists of group
##  numbers in a fixed order, and the data of one of them; see
##  dev/primid-invariants.g, which computes the file.  The data follow from
##  the signs of the groups and the sets <cts> of their cycle types, which
##  PRIMGRP_IdClassData takes from the conjugacy classes.
##
BindGlobal("PRIMGRP_IdClasses", function(deg)
  local PD, keys;
  PD := List([1..NrPrimitiveGroups(deg)], t -> PRIMGrp(deg, t));
  keys := List(PD, e -> [e[2], e[4], Set(e[5]), e[6]]);
  return List(Filtered(Collected(keys), c -> c[2] > 1 and c[1][2] <> "4c"),
              c -> Filtered([1..Length(PD)], i -> keys[i] = c[1]));
end);

BindGlobal("PRIMGRP_IdClassDataFromTypes", function(nrs, signs, cts)
  local all, common, types;

  # a collision between two cycle types of one class could rule out the
  # right group
  all := Union(cts);
  if Length(Set(all, PRIMGRP_CycleTypeHash)) <> Length(all) then
    Error("hash collision for groups ", nrs);
  fi;

  cts := List(cts, ct -> Set(ct, PRIMGRP_CycleTypeHash));
  all := Union(cts);
  common := Intersection(cts);
  types := Difference(all, common);
  return [ nrs, signs, types,
           List(cts, ct -> Sum(Filtered([1..Length(types)],
                                        j -> types[j] in ct), j -> 2^(j-1))) ];
end);

BindGlobal("PRIMGRP_IdClassData", function(deg, nrs)
  local groups;
  groups := List(nrs, i -> PrimitiveGroup(deg, i));
  return PRIMGRP_IdClassDataFromTypes(nrs, List(groups, SignPermGroup),
    List(groups, G -> Set(ConjugacyClasses(G : onlysizes),
                          c -> CycleStructurePerm(Representative(c)))));
end);

#############################################################################
##
#F  PRIMGRP_FilterByCycleTypes( <G>, <PD> )
##
##  <PD> are the library entries that agree with <G> in all fields stored in
##  data/.  Return those not ruled out by the data of data/primid.g (see
##  dev/primid-invariants.g): by the sign of <G>, and by the cycle types of
##  random elements of <G>.  An element whose cycle type a group lacks rules
##  that group out, so the group <G> belongs to is never ruled out, and no
##  library group is built.
##
BindGlobal("PRIMGRP_IDDATA", rec( read := false, data := [] ));

# the data of data/primid.g for degree <deg>, read on first use
BindGlobal("PRIMGRP_IdData", function(deg)
  if not PRIMGRP_IDDATA.read then
    ReadPackage("primgrp", "data/primid.g");
    PRIMGRP_IDDATA.read := true;
  fi;
  if not IsBound(PRIMGRP_IDDATA.data[deg]) then
    return [];
  fi;
  return PRIMGRP_IDDATA.data[deg];
end);

BindGlobal("PRIMGRP_FilterByCycleTypes", function(G, PD)
  local deg, data, pos, s, keep, tries, j, bit;
  deg := NrMovedPoints(G);
  data := First(PRIMGRP_IdData(deg), c -> PD[1][1] in c[1]);
  if data = fail then
    return PD;
  fi;
  pos := List(PD, primdata -> Position(data[1], primdata[1]));
  if fail in pos then
    return PD;
  fi;

  s := SignPermGroup(G);
  keep := Filtered([1..Length(PD)], i -> data[2][pos[i]] = s);

  # while the remaining groups differ in their cycle types
  tries := 0;
  while Length(keep) > 1 and tries < 100
        and Length(Set(keep, i -> data[4][pos[i]])) > 1 do
    tries := tries + 1;
    j := Position(data[3], PRIMGRP_CycleTypeHash(CycleStructurePerm(PseudoRandom(G))));
    if j <> fail then
      bit := 2^(j-1);
      keep := Filtered(keep, i -> IsOddInt(QuoInt(data[4][pos[i]], bit)));
    fi;
  od;
  return PD{keep};
end);

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

  #
  # sign and cycle types of random elements, against data/primid.g
  #
  PD := PRIMGRP_FilterByCycleTypes(grp, PD);
  if Length(PD)=1 then return PD[1][1]; fi;

  #
  # product action: decided outright, on m*k points instead of m^k
  #
  if PD[1][4] = "4c" then
    PD := PRIMGRP_IdentifyProductAction(grp, PD);
    if Length(PD) <> 1 then
      Error("Uh-Oh, this should never happen ", List(PD, primdata -> primdata[1]));
    fi;
    return PD[1][1];
  fi;

  # now we need to create the groups
  cand := List(PD, primdata -> PrimitiveGroup(deg, primdata[1]));

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
