#!/usr/bin/env python3
"""Bring the classical-group names in data/gps*.g to one convention.

The library currently writes the same group several ways: "L(3, 2)" and
"PSL(3, 2)", "PSL(2, 13^2)" and "PSL(2, 169)", "PSU(4,5)" and "PSU(4, 5)".
The majority in each case is PSL-style prefix, a space after the comma, and a
decimal field size, which is also how GAP's own constructors are called.

Only the leading "Xxx(d, q)" of a name is touched; anything after it, such as
":2", ".2_3" or " wreath Sym(2)", is left exactly as it was.

Usage: normalise-names.py [--apply]
"""
import re, sys, glob, collections

PREFIX = {'L': 'PSL'}          # everything else keeps its prefix
HEAD = re.compile(r'^(L|PSL|PGL|PSU|PGU|PSigmaL|PGammaL|PSp|PSO|POmega|Sp|SU|SO)'
                  r'\(\s*(\d+)\s*,\s*([0-9]+(?:\^[0-9]+)?)\s*\)')

def normalise_one(part):
    m = HEAD.match(part.strip())
    if not m:
        return part
    body = part.strip()
    lead = part[:len(part) - len(part.lstrip())]
    trail = part[len(part.rstrip()):]
    pre, d, q = m.group(1), m.group(2), m.group(3)
    pre = PREFIX.get(pre, pre)
    if '^' in q:
        b, e = q.split('^')
        q = str(int(b) ** int(e))
    return lead + f"{pre}({d}, {q})" + body[m.end():] + trail

def normalise(name):
    # a name may list alternatives, "L(2, 2^4):4 = PGammaL(2, 2^4)"
    return '='.join(normalise_one(part) for part in name.split('='))

apply = '--apply' in sys.argv
changes = collections.Counter()
examples = []
for f in sorted(glob.glob('data/gps*.g')):
    s = open(f).read()
    out = []
    pos = 0
    for m in re.finditer(r'"([^"]*)"', s):
        new = normalise(m.group(1))
        if new != m.group(1):
            changes[f] += 1
            if len(examples) < 12:
                examples.append(f"{m.group(1)!r} -> {new!r}")
            out.append(s[pos:m.start()] + '"' + new + '"')
            pos = m.end()
    if out:
        out.append(s[pos:])
        if apply:
            open(f, 'w').write(''.join(out))

total = sum(changes.values())
print(f"{'applied' if apply else 'would change'} {total} names in {len(changes)} files")
for e in examples:
    print("   ", e)
