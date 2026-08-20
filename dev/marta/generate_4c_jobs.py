"""Emit MARTA jobs for the type 4c re-encoding, one per library entry.

    python dev/marta/generate_4c_jobs.py --gaproot '/path/to/gaproot;' \
        --worklist dev/4c-worklist.txt | marta import -

Entries already converted on another machine can be skipped with --skip,
which takes a file of "degree nr" lines.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gaproot", required=True,
                   help="GAP root whose pkg/primgrp is the worktree; keep the ';'")
    p.add_argument("--worklist", required=True,
                   help='file of "degree nr" lines')
    p.add_argument("--skip", help='file of "degree nr" lines to leave out')
    p.add_argument("--max-attempts", type=int, default=2)
    p.add_argument("--count", action="store_true",
                   help="print the number of jobs instead of the jobs")
    args = p.parse_args(argv)

    skip = set()
    if args.skip:
        with open(args.skip) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    skip.add((int(f[0]), int(f[1])))

    n = 0
    with open(args.worklist) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 2:
                continue
            deg, nr = int(f[0]), int(f[1])
            if (deg, nr) in skip:
                continue
            n += 1
            if args.count:
                continue
            job = {
                "case_key": f"primgrp-4c-{deg}-{nr}",
                "command": [
                    "{tools[gap]}", "-q", "-A", "--quitonbreak",
                    "-l", args.gaproot,
                    "-c", f"marta_deg := {deg};; marta_nr := {nr};;",
                    "dev/marta/4c-job.g",
                    "-c", "QUIT;",
                ],
                "params": {"kind": "primgrp-4c", "degree": deg, "nr": nr},
                "slots": 1,
                "tags": ["gap", "primgrp", "4c"],
                "max_attempts": args.max_attempts,
            }
            print(json.dumps(job, sort_keys=True))
    if args.count:
        print(n)
    return 0


if __name__ == "__main__":
    # a downstream `head` closing the pipe is not an error worth a traceback
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        raise SystemExit(0)
