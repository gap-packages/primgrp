#############################################################################
##
##  dev/classical-dual.g
##
##  The type 2 entries that share degree, order, socle and suborbits with a
##  sibling.  Two things can differ and both are searched here:
##
##    - which subgroup of A/T the entry is, when several have order <ext>;
##    - which variety carries the action.  B(2,q) acts on (q^2+1)(q+1)
##      singular points and on just as many totally singular lines, and for
##      odd q those actions are inequivalent.
##
##  Matching is by permutation character -- the multiset over conjugacy
##  classes of (element order, class size, fixed points).  That is an
##  invariant of the permutation group, not of the abstract group, so unlike
##  the class count it separates a group from its dual action.
##
LoadPackage("primgrp");
Read("dev/serialize.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_PermChar", function(G, n)
  return Collected(List(ConjugacyClasses(G),
    c -> [Order(Representative(c)), Size(c), n - NrMovedPoints(Representative(c))]));
end);

##  All candidates at one (ser,d,q,ext): every variety, every subgroup class.
BindGlobal("PRIMGRP_Candidates", function(ser, d, q, ext)
  local out, variety, r, hom, Qg, Q, s, G;
  out := [];
  for variety in ["points", "lines"] do
    r := CALL_WITH_CATCH(PGClassicalNormaliserOn, [ser, d, q, variety]);
    if not r[1] then continue; fi;
    r := r[2];
    if (Size(r.A) / Size(r.T)) mod ext <> 0 then continue; fi;
    if ext = 1 then
      Add(out, rec(variety := variety, w := fail, G := r.T, deg := r.deg));
      continue;
    fi;
    hom := NaturalHomomorphismByNormalSubgroup(r.A, r.T);
    Qg := List(GeneratorsOfGroup(r.A), g -> ImagesRepresentative(hom, g));
    Q := GroupWithGenerators(Qg);
    for s in Filtered(List(ConjugacyClassesSubgroups(Q), Representative),
                      x -> Size(x) = ext) do
      G := PreImage(hom, s);
      Add(out, rec(variety := variety,
                   w := List(GeneratorsOfGroup(s),
                             x -> LetterRepAssocWord(Factorization(Q, x))),
                   G := G, deg := r.deg));
    od;
  od;
  return out;
end);

BindGlobal("PRIMGRP_MatchDual", function(ser, d, q, ext, sibs, out, log)
  local cand, c, nr, stored, inv, hit, i, used, tag;
  cand := PRIMGRP_Candidates(ser, d, q, ext);
  if IsEmpty(cand) then
    for nr in sibs do
      PrintTo(log, ser, "(", d, ",", q, ") ext ", ext, " nr ", nr,
              " SKIP no candidates\n");
    od;
    return;
  fi;
  for c in cand do c.inv := PRIMGRP_PermChar(c.G, c.deg); od;
  used := [];
  for nr in sibs do
    stored := PrimitiveGroup(cand[1].deg, nr);
    inv := PRIMGRP_PermChar(stored, cand[1].deg);
    hit := Filtered([1..Length(cand)], i -> cand[i].inv = inv);
    if Length(hit) <> 1 then
      PrintTo(log, ser, "(", d, ",", q, ") ext ", ext, " nr ", nr,
              " SKIP ", Length(hit), " candidates match\n");
      continue;
    fi;
    if hit[1] in used then
      PrintTo(log, ser, "(", d, ",", q, ") ext ", ext, " nr ", nr,
              " SKIP candidate already taken\n");
      continue;
    fi;
    Add(used, hit[1]);
    c := cand[hit[1]];
    if c.w = fail then
      tag := PRIMGRP_Compact(["cl", ser, d, q, ext, c.variety]);
    else
      tag := PRIMGRP_Compact(["clw", ser, d, q, c.w, c.variety]);
    fi;
    PrintTo(out, c.deg, " ", nr, " ", tag, "\n");
    PrintTo(log, ser, "(", d, ",", q, ") ext ", ext, " nr ", nr,
            " OK ", c.variety, "\n");
  od;
end);
