#############################################################################
##
##  dev/classical-amb.g
##
##  The type 2 natural-action entries whose order does not single them out.
##  Siblings share socle, order, degree and suborbits, so nothing stored tells
##  them apart; the number of conjugacy classes does, and it is an isomorphism
##  invariant, so it can be computed on both sides and matched.
##
##  The match is made here, once.  What gets stored is the word, so reading an
##  entry never has to compute an invariant or search for a subgroup.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
Read("dev/serialize.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_AmbiguousFamily", function(ser, d, q, ext, sibs, out, log)
  local r, hom, Qg, Q, subs, cand, s, G, w, i, stored, ncc, hit, nr, matched;
  r := PGClassicalNormaliser(ser, d, q);
  hom := NaturalHomomorphismByNormalSubgroup(r.A, r.T);
  Qg := List(GeneratorsOfGroup(r.A), g -> ImagesRepresentative(hom, g));
  Q := GroupWithGenerators(Qg);
  subs := Filtered(List(ConjugacyClassesSubgroups(Q), Representative),
                   s -> Size(s) = ext);
  cand := [];
  for s in subs do
    G := PreImage(hom, s);
    Add(cand, rec(ncc := NrConjugacyClasses(G),
                  w := List(GeneratorsOfGroup(s),
                            x -> LetterRepAssocWord(Factorization(Q, x)))));
  od;
  # Match every sibling first.  Two siblings that pick the same candidate are
  # the same abstract group on two different point sets -- points against
  # hyperplanes, say -- which no invariant of the group can separate, and
  # which this encoding cannot express.
  matched := [];
  for nr in sibs do
    stored := PrimitiveGroup(r.deg, nr);
    ncc := NrConjugacyClasses(stored);
    hit := Filtered([1..Length(cand)], i -> cand[i].ncc = ncc);
    if Length(hit) <> 1 then
      Add(matched, fail);
      PrintTo(log, r.deg, " ", nr, " ", ser, "(", d, ",", q, ") ext ", ext,
              " SKIP ", Length(hit), " candidates with ", ncc, " classes\n");
    else
      Add(matched, hit[1]);
    fi;
  od;
  for i in [1..Length(sibs)] do
    if matched[i] = fail then continue; fi;
    if Number(matched, x -> x = matched[i]) > 1 then
      PrintTo(log, r.deg, " ", sibs[i], " ", ser, "(", d, ",", q, ") ext ", ext,
              " SKIP shares its candidate with another entry\n");
      continue;
    fi;
    PrintTo(out, r.deg, " ", sibs[i], " [\"clw\",\"", ser, "\",", d, ",", q, ",",
            PRIMGRP_Compact(cand[matched[i]].w), "]\n");
    PrintTo(log, r.deg, " ", sibs[i], " ", ser, "(", d, ",", q, ") ext ", ext,
            " OK ncc ", cand[matched[i]].ncc, "\n");
  od;
end);
