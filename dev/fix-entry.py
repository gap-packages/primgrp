#!/usr/bin/env python3
"""Replace a token inside one PRIMGRP entry, identified by degree and number.

Needed because the same socle token can be right in one entry and wrong in
another: ["C",[4,3],1] is correct for PrimitiveGroup(3280,1) and wrong for
(3640,1), so a global substitution would fix one and break the other.

Usage: fix-entry.py <deg> <nr> <old> <new> [<deg> <nr> <old> <new> ...]
Refuses to write unless it matched exactly one entry per request.
"""
import sys, glob, re

args = sys.argv[1:]
if not args or len(args) % 4:
    sys.exit(__doc__)

jobs = [tuple(args[i:i+4]) for i in range(0, len(args), 4)]
files = {f: open(f).read() for f in sorted(glob.glob('data/gps*.g'))}

for deg, nr, old, new in jobs:
    hits = 0
    for f, s in files.items():
        m = re.search(r'^PRIMGRP\[%s\]:=' % re.escape(deg), s, re.M)
        if not m:
            continue
        start = m.end()
        end = s.find('\nPRIMGRP[', start)
        end = end if end > 0 else len(s)
        block = s[start:end]
        out = []
        for line in block.split('\n'):
            if re.match(r'\[%s,' % re.escape(nr), line) and old in line:
                line = line.replace(old, new)
                hits += 1
            out.append(line)
        files[f] = s[:start] + '\n'.join(out) + s[end:]
    if hits != 1:
        sys.exit(f"ABORT: degree {deg} nr {nr}: matched {hits} entries, expected 1")
    print(f"  PrimitiveGroup({deg},{nr}): {old} -> {new}")

for f, s in files.items():
    open(f, 'w').write(s)
print(f"{len(jobs)} entries patched")
