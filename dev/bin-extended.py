#!/usr/bin/env python3
"""Bin the staged per-degree files into data/gpsN.g and extend PRIMINDX.

The sweep writes one file per degree; the library wants a few dozen files of
manageable size, indexed by PRIMINDX.  Degrees stay in order and no degree is
split across files, so PRIMINDX[deg] stays a single number.
"""
import sys, os, glob, re

stage, outdir, first_index = sys.argv[1], sys.argv[2], int(sys.argv[3])
limit = int(sys.argv[4]) if len(sys.argv) > 4 else 3_000_000

degs = sorted(int(re.search(r'deg(\d+)\.g$', f).group(1))
              for f in glob.glob(os.path.join(stage, 'deg*.g')))
if not degs:
    sys.exit("no staged degrees found")

index, chunk, size, assign, files = first_index, [], 0, {}, []
def flush():
    global chunk, size, index
    if not chunk:
        return
    name = os.path.join(outdir, f'gps{index}.g')
    with open(name, 'w') as out:
        for d in chunk:
            out.write(open(os.path.join(stage, f'deg{d}.g')).read())
    files.append((index, chunk[0], chunk[-1], os.path.getsize(name)))
    index, chunk, size = index + 1, [], 0

for d in degs:
    n = os.path.getsize(os.path.join(stage, f'deg{d}.g'))
    # a single degree larger than the limit gets its own file
    if chunk and size + n > limit:
        flush()
    assign[d] = index
    chunk.append(d)
    size += n
flush()

total = sum(f[3] for f in files)
print(f"{len(degs)} degrees -> {len(files)} files, {total} bytes")
for i, lo, hi, n in files:
    print(f"  gps{i}.g: {lo}-{hi}  {n/1e6:.2f} MB")

# PRIMINDX entries for the new degrees, formatted like lib/primitiv.grp
lo, hi = degs[0], degs[-1]
vals = [str(assign[d]) for d in range(lo, hi + 1)]
with open(os.path.join(outdir, 'primindx-extension.txt'), 'w') as fh:
    for i in range(0, len(vals), 25):
        fh.write(','.join(vals[i:i+25]) + ',\n')
print(f"PRIMINDX extension for degrees {lo}..{hi} written")
