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
