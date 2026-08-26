#############################################################################
##
##  tst/product4c.g
##
##  Take a group of O'Nan-Scott type 4c apart into wreath product elements.
##
##      gap> ReadPackage("primgrp", "tst/product4c.g");
##      gap> g := PrimitiveGroup(25,24);;
##      gap> r := PRIMGRP_DecomposeProductAction(g);;
##      gap> [ r.m, r.k ];
##      [ 5, 2 ]
##      gap> PRIMGRP_CheckProductAction(g, r);
##      true
##
##  A group of this type acts on the m^k tuples over a set of size m.  Its
##  socle is T^k for a non-abelian simple T, and the orbits of all the factors
##  but the i-th are the fibres of the i-th coordinate; that recovers the
##  identification of the points with tuples, after which each generator reads
##  off as (p_1,...,p_k; sigma).  PGProductAction4c puts one back together.
##
##  This code was instrumental in "compressing" the type 4c entries in data/.
##  We keep it around so we can use it when adding more data, and it also may
##  be helpful for consistency tests in the future.
##
##  Nothing here is read when the package loads.
##

#############################################################################
##
#F  PRIMGRP_OneSocleComponent( <S>, <k> ) . . . . . one factor of a socle T^k
##
##  Return one of the k simple direct factors of <S>, or fail.
##
##  Asking MinimalNormalSubgroups would also do this in theory, but in
##  practice is hopeless once the factors are large: at degree 1089 the socle
##  is A(33)^2 and that call does not finish in twenty minutes, where this
##  function takes seconds.
##
##  A random element of prime order lies with high probability in a socle
##  component. Its S-conjugates then generate that component. Whether it
##  worked can be tested by computing its order, since the k factors are
##  isomorphic and so a component has order |S|^(1/k).
##
##  Thus this is a Las Vegas algorithm. A caller that gets fail should simply
##  ask again.
##
PRIMGRP_OneSocleComponent := function(S, k)
  local x, p, U;
  x := PseudoRandom(S);
  p := Random(PrimeDivisors(Order(x)));
  x := x^(Order(x)/p);
  U := Group(List([1..10], i -> x^PseudoRandom(S)));
  if Size(U)^k <> Size(S) then
    return fail;
  fi;
  return U;
end;

#############################################################################
##
#F  PRIMGRP_SocleComponents( <G>, <S>, <k> ) . . . . all factors of a socle T^k
##
##  Return the k simple direct factors of <S>, the socle of <G>, or fail if
##  the random search above did not settle within its tries.  <G> must act on
##  [ 1 .. LargestMovedPoint(G) ] and permute the factors transitively, which
##  is what type 4c gives.
##
##  <S> normalises each factor, so the others are conjugates of the first
##  under <G> rather than under <S>.  Two are told apart by the partition they
##  cut the points into, which is cheap where comparing the groups themselves
##  is not: the orbits of the i-th factor are the lines in the i-th direction,
##  so distinct factors always differ here.
##
PRIMGRP_SocleComponents := function(G, S, k)
  local deg, U, comps, parts, V, q, tries;

  deg := LargestMovedPoint(G);

  tries := 0;
  repeat
    tries := tries + 1;
    U := PRIMGRP_OneSocleComponent(S, k);
  until U <> fail or tries >= 20;
  if U = fail then
    return fail;
  fi;

  comps := [U];
  parts := [Set(Orbits(U, [1..deg]), Set)];
  tries := 0;
  while Length(comps) < k and tries < 100 do
    tries := tries + 1;
    V := U^PseudoRandom(G);
    q := Set(Orbits(V, [1..deg]), Set);
    if not q in parts then
      Add(comps, V);
      Add(parts, q);
    fi;
  od;
  if Length(comps) <> k then
    return fail;
  fi;

  return comps;
end;

#############################################################################
##
#F  PRIMGRP_DecomposeProductAction( <G>[, <k>] ) . . . . . . . . . . . take apart
##
##  Write <G>, of O'Nan-Scott type 4c on m^k points, over the k coordinates.
##  Returns a record with
##
##      m, k    the size of the inner set and the number of coordinates
##      gens    the generators taken apart, in the order they were taken
##      els     those generators, each as [ p_1, ..., p_k, sigma ]
##      rel     the relabelling of the points that numbers them as tuples
##
##  or a record with a component `err' saying which step did not hold.
##
##  The relabelling is a bijection, so PGProductAction4c(m,k,els) is <G>
##  conjugated by `rel' -- conjugate, not equal.  PRIMGRP_CheckProductAction
##  below tests exactly that, which makes it an identity rather than a
##  comparison of invariants.
##
##  <k> defaults to the width of the socle type, which for an entry of the
##  library is field 8.
##
PRIMGRP_DecomposeProductAction := function(arg)
  local G, k, deg, S, facs, m, i, j, others, orbs, fib, coord, pw, index,
        gens, els, g, sigma, ps, blk, rel, x;

  G := arg[1];
  deg := LargestMovedPoint(G);
  if Length(arg) > 1 then
    k := arg[2];
  else
    k := SocleTypePrimitiveGroup(G).width;
  fi;

  m := RootInt(deg, k);
  if m^k <> deg then
    return rec(err := "degree is not a k-th power");
  fi;

  S := Socle(G);
  facs := PRIMGRP_SocleComponents(G, S, k);
  if facs = fail then
    return rec(err := "could not split the socle into k factors");
  fi;

  # the fibres of coordinate i are the orbits of everything but the i-th factor
  fib := [];
  for i in [1..k] do
    others := Subgroup(S, Concatenation(List(Filtered([1..k], j -> j <> i),
                                             j -> GeneratorsOfGroup(facs[j]))));
    orbs := Set(List(Orbits(others, [1..deg]), Set));
    if Length(orbs) <> m then
      return rec(err := "wrong number of fibres");
    fi;
    fib[i] := orbs;
  od;

  # the point relabelling, c_1 most significant
  coord := List([1..deg], x -> List([1..k],
             i -> First([1..m], a -> x in fib[i][a])));
  pw := List([1..k], i -> m^(k-i));
  index := c -> 1 + Sum([1..k], i -> (c[i]-1)*pw[i]);
  rel := PermList(List([1..deg], x -> index(coord[x])));

  # each generator permutes the fibres, giving sigma, and acts on each
  # coordinate, giving p_i
  gens := GeneratorsOfGroup(G);
  if Length(gens) > 2 then
    gens := SmallGeneratingSet(G);      # dear, and cannot beat two
  fi;
  els := [];
  for g in gens do
    sigma := [];
    ps := [];
    for i in [1..k] do
      blk := Set(OnTuples(fib[i][1], g));
      j := First([1..k], j -> ForAny(fib[j], f -> f = blk));
      if j = fail then
        return rec(err := "a generator does not permute the fibres");
      fi;
      sigma[i] := j;
      ps[i] := PermList(List([1..m],
                 a -> First([1..m], b -> fib[j][b] = Set(OnTuples(fib[i][a], g)))));
    od;
    Add(els, Concatenation(ps, [PermList(sigma)]));
  od;

  return rec(m := m, k := k, els := els, rel := rel, gens := gens);
end;

#############################################################################
##
#F  PRIMGRP_CheckProductAction( <G>, <r> ) . . . . . . . . put back together
##
##  Is the group PGProductAction4c builds from the decomposition <r> the group
##  <G> relabelled?  An identity, so a true here needs no further checking.
##
##  The generators are compared one by one, in the order they were taken apart,
##  rather than the two groups being compared as groups.  It says more -- each
##  element rebuilds, not merely the group they generate -- and it says it
##  without a stabiliser chain, which at these degrees and orders is the
##  difference between seconds and hours.
##
PRIMGRP_CheckProductAction := function(G, r)
  if IsBound(r.err) then
    return false;
  fi;
  return GeneratorsOfGroup(PGProductAction4c(r.m, r.k, r.els))
         = List(r.gens, x -> x^r.rel);
end;
