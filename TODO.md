# TODO

## Documentation of the degree range

`doc/manual.xml:58` still says the library covers degree &lt; 4096. Since
degrees 4096 to 8191 moved into `data/`, that is wrong however the package is
eventually shipped -- but what it should say depends on a decision not yet
made: which data files go in the release tarball and which are fetched as
ArtifactManager artifacts. Write the sentence once that is settled, so the
manual does not have to be corrected twice.

`README.md` has the same sentence and is already correct (degree &lt; 8192);
it no longer mentions a separate download.

## ArtifactManager

Splitting `data/` between the tarball and downloadable artifacts is separate
work, to be done alongside the ArtifactManager integration. Nothing in the
package refers to artifacts at the moment; `.release` gzips all of `data/*.g`.

## How to write the field size in a name

Names now say `PSL(2, 169)`, not `PSL(2, 13^2)`. That was chosen because it is
the majority convention in the existing data (530 against 200), it matches how
GAP's own constructors are called, and having one convention was worth more
than either particular choice.

But the prime-power form carries information the decimal does not: whether q
is prime, a square, a higher power. That governs which extensions exist --
PSigmaL is only larger than PSL when q is not prime -- so it is exactly what a
reader wants to see. Worth revisiting.

The right moment is once every name comes from a constructor rather than from
the data files. Then the convention lives in one place, `PRIMGRP_FieldString`
in `lib/primitiv.gi` (removed for now, but that is where it belongs), and the
choice is a one-line change instead of a 230-entry rewrite. As of now the
classical families are constructor-generated, but the sporadic, alternating
and affine names are not.

## What is still stored as permutations, and how it might be encoded

At the time of writing 1121 entries hold explicit permutations, 32.7 MB of
`data/`, in five groups.  Each is listed with what is known about it, so the
next attempt does not start from the same measurements.

### Type 4c, 142 entries, 9.76 MB

Outstanding MARTA jobs at the time of the freeze.  The machinery works --
1412 of 1553 are already converted -- so these need running, not inventing.
One, degree 5041 entry 208, failed twice on tutulla and wants a look; one,
degree 2401 entry 1173, hung for hours locally and is excluded by the skip
list.

### Type 2 on cosets of a maximal subgroup, 950 entries, 21.9 MB

Not "non-natural actions" -- that label was wrong.  For socle PSL(2,q) the
maximal subgroups are classified by Dickson, and all 232 such entries sit on
one of: Borel (the natural action), dihedral D_2(q+1)/d or D_2(q-1)/d, A4, S4,
A5, or a subfield subgroup.  The other socle families are equally classified.

The obstacle is cost, not theory.  Building the group as a coset action was
measured at degree 8128:

    MaximalSubgroupClassReps + FactorCosetAction   47 s
    direct normaliser + FactorCosetAction         170 ms
    reading the stored permutations                 7 ms

and of that 170 ms, 154 is FactorCosetAction itself, so writing the normaliser
generators down universally -- worth doing anyway, since the random search for
a torus element is not reproducible -- only saves 23 ms.  Multiplied over a
thousand entries the difference between 7 s and 170 s of library time is
noticeable.

The way through is a direct construction for the actions that are geometric or
combinatorial, computing each point's image arithmetically instead of looking
it up, exactly as `PRIMGRP_AffineAction` does:

  - `D_2(q-1)/d` has degree Binomial(q+1,2) and *is* the action on 2-subsets
    of the projective line, so `PGOnSets` already covers it.
  - `D_2(q+1)/d`, 102 entries, 2.54 MB, has degree q(q-1)/2 counting the
    Frobenius-conjugate pairs of non-rational points of PG(1,q^2).  Clean, and
    the largest single family.
  - subfield subgroups, 28 entries, 0.67 MB, act on Baer sublines.
  - A4, S4 and A5, 36 entries, 0.78 MB, have no such description.  A coset
    action is the honest encoding and these should probably stay as they are.

### Types 4a and 4b, 24 entries, 0.91 MB

Two diagonal components in product action: socle T^4 of degree |T|^2, so
3600 = 60^2 for T = Alt(5).  `PRIMGRP_Extract4c` assumes the socle's k factors
give m^k = degree and so cannot see this; it needs generalising to components
that are themselves diagonal rather than simple.

### Type 3b of socle width 3, 5 entries, 0.17 MB

Degree |T|^2, the action being on cosets of the diagonal in T^3 rather than on
the elements of T.  `PGDiagonal` handles width 2 only.

### Classical entries the matcher has not settled

Around 70 entries sit at their socle's natural action degree but are still
permutations, because the ambiguity matcher has not reached them -- twice now
the family list was sorted so that the largest group, socle L, came last and
the run was stopped before it got there.  Order the work by payoff.

Genuinely refused, rather than unreached:

  - `D(4,3)` and `2D(4,3)`, where triality gives PSO+(8,3) three inequivalent
    actions on 1120 points and the normaliser built from the singular points
    is too small by half;
  - `PSL(5,2).2`, `PSL(6,2):2`, `PSL(7,2).2` on 2-subsets, where the graph
    automorphism maps points to hyperplanes and so is not realised on points.
    A "hyperplanes" variety in `PGClassicalNormaliserOn`, alongside the
    existing "points" and "lines", would settle these.

## Performance is a property to test, not to assume

Every conversion in this work was checked for correctness and none for cost,
and `PrimitiveGroup(6561,1)` reached the maintainer taking long enough to
interrupt.  The cause was a linear scan over a hand-built point list;
`tst/testinstall/performance.tst` now guards it.  Two things follow:

  - a construction-time sweep over all 8190 degrees would catch the next one,
    and fits the MARTA job shape already in `dev/marta/`;
  - the cost of `AllPrimitiveGroups` and `PrimitiveGroupsIterator` over a wide
    range has never been measured, and the lazy entries and constructors are
    exactly the kind of change that could hurt it.

## Make PrimitiveIdentification cheaper by storing invariants

`PrimitiveIdentification` does not fall back on `IsomorphismGroups`, contrary
to a reasonable guess. It filters on a ladder of invariants: size, O'Nan-Scott
type, suborbits and transitivity come free from the data file, and then, past
the comment "now we need to create the groups", it computes the socle-quotient
`IdGroup`, `AbelianInvariants`, Sylow subgroup orbit shapes, lower central
series of those Sylow subgroups, Frattini subgroup orders, and finally the
cycle structures of conjugacy class representatives.

The cost is not the invariants. It is that everything after the fourth is
computed on *constructed* candidates: to identify one group of degree 7776 it
builds a dozen groups of degree 7776 and takes Sylow subgroups and conjugacy
classes of each. Adding more invariants would make this worse, not better.

The fix is to precompute a few of them once and store them next to the entry,
so the filter runs on stored data for longer and constructs candidates only as
a last resort. `AbelianInvariants` and the socle-quotient `IdGroup` are the
obvious first two: small, cheap to store, and they cut deeply. A separate
table, as with `PRIMGRP_SIMSNO`, avoids touching the entry format.

This is what blocks compressing the product-action groups above degree 4095 --
38.9 MB, the largest single block left -- because each conversion has to be
confirmed by an identification, and at those degrees each one costs minutes.
