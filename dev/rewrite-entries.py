#!/usr/bin/env python3
"""Replace whole entries in data/gps*.g with the text given for each.

Input is lines of "<degree> <nr> <text>".  The entry's whole bracketed span is
replaced by <text>, so an entry becomes the call that produces it.

Companion to dev/apply-entries.py, which replaces field 9 and leaves the entry
in place; this replaces the entry itself.  Both refuse to write unless every
request matched exactly one entry, and both record the degrees they touched.

The files are not one entry to a line -- a degree's list wraps and its entries
are indented into it -- so this walks brackets over the file text.

Usage: rewrite-entries.py <datadir> <worklist>
"""
import sys, os, glob, re

DIRTY = 'dev/sweep-dirty.txt'


def scan_top(s, start):
    """(begin, end) of each top-level [...] inside the list opening at start"""
    assert s[start] == '['
    i, depth, quoted, begin = start + 1, 0, False, None
    while i < len(s):
        ch = s[i]
        if quoted:
            if ch == '\\':
                i += 2
                continue
            if ch == '"':
                quoted = False
        elif ch == '"':
            quoted = True
        elif ch == '[':
            if depth == 0:
                begin = i
            depth += 1
        elif ch == ']':
            depth -= 1
            if depth == 0:
                yield (begin, i + 1)
            elif depth < 0:
                return
        i += 1


def first_field(entry):
    """the entry number, read off the text"""
    m = re.match(r'\[\s*(\d+)\s*,', entry)
    return int(m.group(1)) if m else None


def main(argv):
    if len(argv) != 2:
        sys.exit(__doc__)
    datadir, worklist = argv

    want = {}
    for line in open(worklist):
        line = line.rstrip('\n')
        if not line.strip() or line.startswith('#'):
            continue
        deg, nr, text = line.split(' ', 2)
        key = (int(deg), int(nr))
        if key in want:
            sys.exit(f"ABORT: {deg} {nr} requested twice")
        want[key] = text

    hits, edits = {}, {}
    for path in sorted(glob.glob(os.path.join(datadir, 'gps*.g'))):
        s = open(path).read()
        cuts = []
        for m in re.finditer(r'^PRIMGRP\[(\d+)\]:=\s*', s, re.M):
            deg = int(m.group(1))
            if not any(d == deg for d, _ in want):
                continue
            lb = s.index('[', m.end() - 1)
            for (b, e) in scan_top(s, lb):
                nr = first_field(s[b:e])
                if nr is None or (deg, nr) not in want:
                    continue
                cuts.append((b, e, want[(deg, nr)]))
                hits[(deg, nr)] = hits.get((deg, nr), 0) + 1
        if cuts:
            out, last = [], 0
            for (b, e, text) in cuts:
                out.append(s[last:b])
                out.append(text)
                last = e
            out.append(s[last:])
            edits[path] = ''.join(out)

    missing = [k for k in want if hits.get(k, 0) != 1]
    if missing:
        for deg, nr in sorted(missing)[:5]:
            print(f"  degree {deg} nr {nr}: matched {hits.get((deg, nr), 0)},"
                  " expected 1", file=sys.stderr)
        sys.exit(f"ABORT: {len(missing)} requests did not match exactly one entry")

    for path, text in edits.items():
        open(path, 'w').write(text)

    touched = sorted({deg for deg, _ in want})
    seen = set()
    if os.path.exists(DIRTY):
        seen = {int(x) for x in open(DIRTY).read().split()}
    with open(DIRTY, 'a') as fh:
        for d in touched:
            if d not in seen:
                fh.write(f"{d}\n")
    print(f"{len(want)} entries rewritten across {len(touched)} degrees "
          f"in {len(edits)} files")


main(sys.argv[1:])
