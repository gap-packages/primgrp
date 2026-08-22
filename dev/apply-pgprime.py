"""Replace the affine entries at prime degrees with PGPrime(<d>).

For prime p the affine primitive groups are the subgroups of AGL(1,p)
containing the translations, one per divisor of p-1, so p and d determine the
entry completely.  The block is contiguous but does not always start at nr 1:
at degree 263 it is PGAlt, PGSym, then nrs 3..6.

Entry lines are located by bracket balance.  Stripping trailing ']' instead
mangles the last entry of a degree, which ends ']];'.
"""
import re, glob, sys

def entry_span(line):
    if not line.startswith('['):
        return None
    d = q = i = 0
    q = False
    while i < len(line):
        ch = line[i]
        if q:
            if ch == '\\': i += 2; continue
            if ch == '"': q = False
        elif ch == '"': q = True
        elif ch in '[(': d += 1
        elif ch in ')]':
            d -= 1
            if d == 0: return i + 1
        i += 1
    return None

def split_top(s):
    out, d, cur, q = [], 0, [], False
    i = 0
    while i < len(s):
        ch = s[i]
        if q:
            cur.append(ch)
            if ch == '\\': cur.append(s[i+1]); i += 2; continue
            if ch == '"': q = False
            i += 1; continue
        if ch == '"': q = True; cur.append(ch)
        elif ch in '[(': d += 1; cur.append(ch)
        elif ch in '])': d -= 1; cur.append(ch)
        elif ch == ',' and d == 0: out.append(''.join(cur)); cur = []
        else: cur.append(ch)
        i += 1
    out.append(''.join(cur))
    return out

def isprime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0: return False
        f += 2
    return True

def divisors(n):
    small, large = [], []
    i = 1
    while i * i <= n:
        if n % i == 0:
            small.append(i)
            if i != n // i: large.append(n // i)
        i += 1
    return small + large[::-1]

def main():
    files = sorted(glob.glob('data/gps*.g'),
                   key=lambda s: int(re.search(r'\d+', s).group()))
    replaced = saved = 0
    degrees = []
    for path in files:
        lines = open(path).read().split('\n')
        deg = None
        block = []            # (index, nr) of affine entries of the current degree
        out = list(lines)
        changed = False

        def flush(deg, block):
            nonlocal replaced, saved, changed
            if deg is None or not isprime(deg) or deg < 5:
                return
            D = divisors(deg - 1)
            if len(block) != len(D):
                return
            nrs = [b[1] for b in block]
            if nrs != list(range(nrs[0], nrs[0] + len(D))):
                return
            for (idx, _), d in zip(block, D):
                line = lines[idx]
                e = entry_span(line)
                tail = line[e:]
                saved += len(line) - (len('PGPrime(%d)' % d) + len(tail))
                out[idx] = 'PGPrime(%d)%s' % (d, tail)
                replaced += 1
            degrees.append(deg)
            changed = True

        for i, line in enumerate(lines):
            m = re.match(r'PRIMGRP\[(\d+)\]', line)
            if m:
                flush(deg, block)
                deg = int(m.group(1)); block = []
                continue
            if deg is None or not line.startswith('['):
                continue
            e = entry_span(line)
            if e is None: continue
            fs = split_top(line[1:e-1])
            if len(fs) < 9: continue
            try: nr = int(fs[0])
            except ValueError: continue
            if fs[3] == '"1"' and fs[7] == '["Z",%d,1]' % deg:
                block.append((i, nr))
        flush(deg, block)
        if changed:
            open(path, 'w').write('\n'.join(out))
    print("entries replaced: %d over %d prime degrees" % (replaced, len(degrees)))
    print("bytes saved: %d" % saved)
    with open('dev/pgprime-degrees.txt', 'w') as fh:
        for d in sorted(degrees): fh.write("%d\n" % d)

if __name__ == '__main__':
    sys.exit(main())
