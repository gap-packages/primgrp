[![CI](https://github.com/gap-packages/primgrp/actions/workflows/CI.yml/badge.svg)](https://github.com/gap-packages/primgrp/actions/workflows/CI.yml)
[![Code Coverage](https://codecov.io/github/gap-packages/primgrp/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/primgrp)

# GAP Primitive Groups Library

The PrimGrp package provides the library of primitive permutation 
groups which includes, up to permutation isomorphism (i.e., up to
conjugacy in the corresponding symmetric group), all primitive 
permutation groups of degree < 8192. Groups of degree 4096 to 8191 are about
1 GB and are not shipped with the package; they have to be obtained
separately, in one of two ways.

If the [ArtifactManager](https://github.com/gap-packages/ArtifactManager)
package is available, it will do it for you:

```gap
gap> LoadPackage("ArtifactManager");
gap> FetchArtifact("primgrp", "degree4096to8191");
```

It downloads the data, checks it against a SHA256 checksum, and caches it
outside the package directory — so the data survives reinstalling PrimGrp, and
`RemoveArtifact("primgrp", "degree4096to8191")` reclaims the space again. The
download is not started automatically, because 1 GB is not something a call to
`PrimitiveGroup` should decide on your behalf.

Otherwise, download the data by hand from
<https://doi.org/10.5281/zenodo.10411366> and follow the instructions given
there. Data installed that way is still found, and takes precedence.

It has been formerly a part of the core GAP system, and has been
converted to a separate GAP package in 2017 for the GAP 4.9 release.

## Documentation

Full information and documentation can be found in the manual, available
as PDF `doc/manual.pdf` or as HTML `doc/chap0.html`, or on the package
homepage at

  <https://gap-packages.github.io/primgrp/>


## Bug reports and feature requests

Please submit bug reports and feature requests via our GitHub issue tracker:

  <https://github.com/gap-packages/primgrp/issues>


# License

PrimGrp is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

For details see the LICENSE file.

