#############################################################################
##
##  dev/product4c.g
##
##  Rewrite a type 4c entry as elements of Sym(m) wreath Sym(k).
##
##  The socle is T^k, and the orbits of all factors but the i-th are the
##  fibres of the i-th coordinate.  That recovers the identification of the
##  m^k points with k-tuples, after which each generator reads off as
##  (p_1,...,p_k; sigma).
##
##  The relabelling of points is a bijection, so the group built back from
##  the tuples is the original conjugated by it -- conjugate, not equal, which
##  is what this phase allows.  The check below is on that conjugate, so it is
##  an identity rather than an invariant comparison.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
Read("dev/serialize.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_Extract4c", function(G, deg, k)
  local S, facs, m, i, others, orbs, fib, coord, pw, index, gens, els, g,
        sigma, ps, a, b, blk, j, rel, x, c;
  S := Socle(G);
  facs := MinimalNormalSubgroups(S);
  if Length(facs) <> k then return rec(err := "socle width"); fi;
  m := RootInt(deg, k);
  if m^k <> deg then return rec(err := "degree not a k-th power"); fi;
  fib := [];
  for i in [1..k] do
    others := Subgroup(S, Concatenation(List(Filtered([1..k], j -> j <> i),
                                             j -> GeneratorsOfGroup(facs[j]))));
    orbs := Set(List(Orbits(others, [1..deg]), Set));
    if Length(orbs) <> m then return rec(err := "fibre count"); fi;
    fib[i] := orbs;
  od;
  coord := List([1..deg], x -> List([1..k],
             i -> First([1..m], a -> x in fib[i][a])));
  pw := List([1..k], i -> m^(k-i));
  index := c -> 1 + Sum([1..k], i -> (c[i]-1)*pw[i]);
  rel := PermList(List([1..deg], x -> index(coord[x])));
  gens := SmallGeneratingSet(G);
  els := [];
  for g in gens do
    sigma := []; ps := [];
    for i in [1..k] do
      blk := Set(OnTuples(fib[i][1], g));
      j := First([1..k], j -> ForAny(fib[j], f -> f = blk));
      if j = fail then return rec(err := "generator does not permute fibres"); fi;
      sigma[i] := j;
      ps[i] := PermList(List([1..m],
                 a -> First([1..m], b -> fib[j][b] = Set(OnTuples(fib[i][a], g)))));
    od;
    Add(els, Concatenation(ps, [PermList(sigma)]));
  od;
  return rec(m := m, k := k, els := els, rel := rel, gens := gens);
end);

##  Convert one entry, and insist the rebuild is the original relabelled.
BindGlobal("PRIMGRP_Convert4c", function(deg, nr)
  local G, e, k, r, H;
  e := PRIMGrp(deg, nr);
  k := e[8][3];
  G := PrimitiveGroup(deg, nr);
  r := PRIMGRP_Extract4c(G, deg, k);
  if IsBound(r.err) then return rec(err := r.err); fi;
  H := PGProductAction4c(r.m, r.k, r.els);
  if H <> Image(ConjugatorIsomorphism(G, r.rel), G) then
    return rec(err := "rebuild differs from the relabelled original");
  fi;
  return rec(text := PRIMGRP_Compact(["4c", r.m, r.k, r.els]));
end);
