#############################################################################
##
##  dev/diagonal.g
##
##  O'Nan-Scott 3a and 3b of socle width 2: T x T on the |T| elements of T.
##  Everything between T^2 and N = T^2.(Out(T) x 2) is such an entry, and
##  |N/T^2| is 2 to 8, so the choice is among a handful of subgroups.
##
##  Several entries of one degree can share an order -- degree 60 has three
##  subgroups of order 2, one of type 3a and two of type 3b -- so they are
##  told apart by permutation character, which is an invariant of the action.
##
LoadPackage("primgrp");
ReadPackage("primgrp", "tst/testutils.g");
Read("dev/serialize.g");
SetInfoLevel(InfoWarning, 0);

BindGlobal("PRIMGRP_DiagPermChar", function(G, n)
  return Collected(List(ConjugacyClasses(G),
    c -> [Order(Representative(c)), Size(c), n - NrMovedPoints(Representative(c))]));
end);

##  All entries of one degree at once: one normaliser, one subgroup lattice.
BindGlobal("PRIMGRP_MatchDiagonal", function(tspec, sibs, out, log)
  local r, hom, Qg, Q, cand, s, G, ext, nr, stored, inv, hit, used, i, tag;
  r := PGDiagonalNormaliser(tspec);
  hom := NaturalHomomorphismByNormalSubgroup(r.N, r.T2);
  Qg := List(GeneratorsOfGroup(r.N), g -> ImagesRepresentative(hom, g));
  Q := GroupWithGenerators(Qg);
  cand := [];
  for s in List(ConjugacyClassesSubgroups(Q), Representative) do
    G := PreImage(hom, s);
    Add(cand, rec(ext := Size(s),
                  w := List(GeneratorsOfGroup(s),
                            x -> LetterRepAssocWord(Factorization(Q, x))),
                  nclass := Number(List(ConjugacyClassesSubgroups(Q), Representative),
                                   x -> Size(x) = Size(s)),
                  inv := PRIMGRP_DiagPermChar(G, r.deg)));
  od;
  used := [];
  for nr in sibs do
    stored := PrimitiveGroup(r.deg, nr);
    inv := PRIMGRP_DiagPermChar(stored, r.deg);
    hit := Filtered([1..Length(cand)], i -> cand[i].inv = inv);
    if Length(hit) <> 1 then
      PrintTo(log, r.deg, " ", nr, " SKIP ", Length(hit), " candidates match\n");
    elif hit[1] in used then
      PrintTo(log, r.deg, " ", nr, " SKIP candidate already taken\n");
    else
      Add(used, hit[1]);
      if cand[hit[1]].nclass = 1 then
        tag := PRIMGRP_Compact(["diag", tspec, cand[hit[1]].ext]);
      else
        tag := PRIMGRP_Compact(["diagw", tspec, cand[hit[1]].w]);
      fi;
      PrintTo(out, r.deg, " ", nr, " ", tag, "\n");
      PrintTo(log, r.deg, " ", nr, " OK ext ", cand[hit[1]].ext,
              " classes ", cand[hit[1]].nclass, "\n");
    fi;
  od;
end);
