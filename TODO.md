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
