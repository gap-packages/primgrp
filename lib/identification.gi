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
#F  PRIMGRP_AffinePointStabilizer( <G> )
##
##  For a primitive group <G> of affine type, a record describing a point
##  stabiliser: a faithful permutation representation of small degree, and
##  what is needed to write its elements as matrices over the socle.
##
BindGlobal("PRIMGRP_AffinePointStabilizer", function(G)
  local dom, stab, orbs, hom;
  dom := MovedPoints(G);
  stab := Stabilizer(G, dom[1]);

  # the stabiliser acts irreducibly on the socle, so an orbit of non-zero
  # vectors spans it, and the action on any such orbit is faithful
  orbs := Orbits(stab, dom{[2..Length(dom)]});
  hom := ActionHomomorphism(stab, orbs[PositionMinimum(orbs, Length)],
                            "surjective");

  return rec( stab := stab, hom := hom, perm := Image(hom),
              pcgs := Pcgs(Socle(G)),
              field := GF(Factors(Length(dom))[1]) );
end);

#############################################################################
##
#F  PRIMGRP_IsConjugateAffine( <a>, <b> )
##
##  Two primitive groups of affine type and equal degree are conjugate in the
##  symmetric group if and only if their point stabilisers are conjugate in
##  GL(d,p): that is, if some isomorphism between the stabilisers carries the
##  natural module of one to that of the other.  An inner automorphism does not
##  change a module up to isomorphism, so one isomorphism for each class of
##  outer automorphisms is tried.
##
##  Searching GL(d,p) for a conjugating element instead can run for hours
##  when there is none.
##
##  <a> and <b> are the records PRIMGRP_AffinePointStabilizer returns for the
##  two groups.  They are arguments so that a caller comparing one group with
##  many computes its record, and the automorphism group stored in a.perm,
##  once.
##
BindGlobal("PRIMGRP_IsConjugateAffine", function(a, b)
  local iso, gens, module, m, aut, inn;

  # cyclic groups of prime order
  if IsTrivial(a.perm) then
    return IsTrivial(b.perm);
  fi;

  iso := IsomorphismGroups(a.perm, b.perm);
  if iso = fail then
    return false;
  fi;

  # the natural module of the stabiliser in <r>, with the generators in <elms>
  module := function(r, elms)
    elms := List(elms, x -> PreImagesRepresentative(r.hom, x));
    return GModuleByMats(LinearActionLayer(r.stab, elms, r.pcgs), r.field);
  end;

  gens := GeneratorsOfGroup(a.perm);
  m := module(a, gens);
  aut := AutomorphismGroup(a.perm);
  inn := InnerAutomorphismsAutomorphismGroup(aut);
  return ForAny(RightTransversal(aut, inn), alpha ->
    MTX.IsomorphismModules(m,
      module(b, List(gens, x -> Image(iso, Image(alpha, x))))) <> fail);
end);

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
  x := PseudoRandom(G)^e;
  if IsOne(x) then
    return fail;
  fi;
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
#F  PRIMGRP_ProductActionSize( <G>, <PD> )
##
##  The order of <G>, computed in its action on the m*k pairs, if <G> is of
##  type 4c on m^k points; otherwise, or if the random search in
##  PRIMGRP_WreathElements misses, fail.  <PD> holds the library entries of the
##  degree of <G>.  Only socles Alt(m)^k are tried, the case where Size on m^k
##  points is slow.
##
##  The action on the pairs is faithful, and PRIMGRP_WreathElements verifies
##  that it is the action of <G>, so the order is exact.
##
BindGlobal("PRIMGRP_ProductActionSize", function(G, PD)
  local deg, k, m, e, i, els;
  deg := NrMovedPoints(G);
  PD := Filtered(PD, primdata -> primdata[4] = "4c" and primdata[8][1] = "A"
                       and primdata[8][2]^primdata[8][3] = deg);
  for k in Set(PD, primdata -> primdata[8][3]) do
    m := RootInt(deg, k);

    # a multiple of the index of the socle Alt(m)^k
    e := Lcm(List(Filtered(PD, primdata -> primdata[8][3] = k),
                  primdata -> primdata[2] / (Factorial(m)/2)^k));

    # a try costs 10 ms and succeeds in a quarter to a half of the cases for
    # Alt(m)^2, m = 68..89; Size on the m^k points takes minutes
    for i in [1..50] do
      els := PRIMGRP_WreathElements(G, m, k, e);
      if els <> fail then
        return Size(Group(List(els, x -> PRIMGRP_WreathToPerm(m, k, x))));
      fi;
    od;
  od;
  return fail;
end);

#############################################################################
##
#F  PRIMGRP_StarAction( <G>, <PD> )
##
##  The action of <G> on m stars, if <G> acts like a group with socle Alt(m)
##  on the 2-subsets of [ 1 .. m ]; otherwise, or if the random search misses,
##  fail.  <PD> holds the library entries of the degree of <G>.  Size on the
##  2-subsets takes 50 s for Sym(111), against 15 ms on the stars.
##
##  The star of i is the set of the m-1 pairs containing i.  If <G> permutes m
##  sets so that every point lies in exactly two and no two points in the same
##  two, then a point is the intersection of its two sets, so <G> is its
##  action on the sets, acting on their 2-subsets.
##
BindGlobal("PRIMGRP_StarAction", function(G, PD)
  local dom, deg, m, gens, a, sv, queue, j, i, x, rep, H, delta, star, O, in2,
        try;
  dom := MovedPoints(G);
  deg := Length(dom);
  PD := First(PD, primdata -> primdata[4] = "2" and primdata[8][1] = "A"
                     and Binomial(primdata[8][2], 2) = deg);
  if PD = fail then
    return fail;
  fi;
  m := PD[8][2];

  # Schreier vector of the orbit of a: sv[x] is the generator that reached x
  gens := GeneratorsOfGroup(G);
  a := dom[1];
  sv := [];
  sv[a] := 0;
  queue := [a];
  for j in queue do
    for i in [1..Length(gens)] do
      x := j^gens[i];
      if not IsBound(sv[x]) then
        sv[x] := i;
        Add(queue, x);
      fi;
    od;
  od;

  # an element mapping a to x
  rep := function(x)
    local u;
    u := ();
    while x <> a do
      u := gens[sv[x]] * u;
      x := x / gens[sv[x]];
    od;
    return u;
  end;

  for try in [1..10] do
    # random elements of the stabiliser of a = {i, j}, which has the orbit
    # delta of pairs that meet {i, j} in one point; the pairs containing i
    # form a block of size m-2 in it, found from delta[1] and any other point
    # of it, one in two
    H := Group(List([1..10], k -> PseudoRandom(G)), ());
    H := Group(List(GeneratorsOfGroup(H), g -> g / rep(a^g)), ());
    delta := First(Orbits(H, dom), o -> Length(o) = 2*(m-2));
    if delta = fail then
      continue;
    fi;
    star := fail;
    for x in delta{[2..Length(delta)]} do
      star := Blocks(H, delta, [delta[1], x])[1];
      if Length(star) = m-2 then
        break;
      fi;
      star := fail;
    od;
    if star = fail then
      continue;
    fi;

    O := Orbit(G, Set(Concatenation([a], star)), OnSets);
    if Length(O) <> m then
      continue;
    fi;
    in2 := [];
    for x in dom do
      in2[x] := [];
    od;
    for i in [1..m] do
      for x in O[i] do
        Add(in2[x], i);
      od;
    od;
    in2 := in2{dom};
    if ForAll(in2, l -> Length(l) = 2) and Length(Set(in2)) = deg then
      # Size alone takes 11 s for Alt(128); the giant tests set it
      H := Action(G, O, OnSets);
      if not IsNaturalAlternatingGroup(H) then
        IsNaturalSymmetricGroup(H);
      fi;
      return H;
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
local dom,deg,PD,s,cand,a,p,f,q,gl;
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

  # get information of all primitive groups of the given degree
  PD := List([1 .. NrPrimitiveGroups(deg)], t -> PRIMGrp(deg, t));

  # Size by Schreier-Sims takes minutes for some groups of large degree: 143 s
  # for D(2*8191), inside the tests for Alt(n) and Sym(n) below, and 60 to
  # 120 s for conjugates of 5329 #280, with socle Alt(73)^2.  A solvable group
  # gets its order from the pcgs IsSolvableGroup finds, which costs at most
  # 55 ms on the others, a product action from its action on m*k points, and
  # a 2-subset action from that on m points.
  if not HasSize(grp) and not IsSolvableGroup(grp) then
    s := PRIMGRP_ProductActionSize(grp, PD);
    if s = fail then
      a := PRIMGRP_StarAction(grp, PD);
      if a <> fail then
        s := Size(a);

        # grp is Alt(m) or Sym(m) on 2-subsets, which is primitive for m >= 5
        # and determined by its order; IsPrimitive takes up to 11 s for
        # Alt(114) on 6441 points
        if IsNaturalAlternatingGroup(a) or IsNaturalSymmetricGroup(a) then
          return First(PD, primdata -> primdata[4] = "2"
                       and primdata[8][1] = "A"
                       and primdata[8][2] = NrMovedPoints(a)
                       and primdata[2] = s)[1];
        fi;
      fi;
    fi;
    if s <> fail then
      SetSize(grp, s);
    fi;
  fi;

  if not HasSize(grp) then
    if IsNaturalAlternatingGroup(grp) then
      SetSize(grp, Factorial(deg)/2);
    elif IsNaturalSymmetricGroup(grp) then
      SetSize(grp, Factorial(deg));
    fi;
  fi;

  if not (IsTransitive(grp,dom) and IsPrimitive(grp,dom)) then
    Error("Group must operate primitively");
  fi;

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
    a := PRIMGRP_AffinePointStabilizer(grp);

    # grp is conjugate to one of cand, so the last needs no test
    p := PositionProperty(cand{[1..Length(cand)-1]}, g ->
      PRIMGRP_IsConjugateAffine(a, PRIMGRP_AffinePointStabilizer(g)));
    if p = fail then
      p := Length(cand);
    fi;
    cand := [ cand[p] ];
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
