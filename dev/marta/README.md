# Running the long primgrp jobs under MARTA

Two jobs on this branch outgrew a laptop: the consistency sweep over degrees
2..8191, and the type 4c re-encoding of 1553 entries.  Both are one-unit-per-
job and independent, so they suit MARTA directly.

Local cost, one core: the sweep reached 4869 of 8190 degrees in 20.5 h, degree
4096 alone taking 2.55 h; 4c did 729 of 1553 entries in 13 h, of which 7.2 h
was a single entry that hung.  Cost per unit is very uneven, which is the
argument for one job per unit rather than per block.

## Worker requirements

- `--tool gap=<path>`, and `transgrp` loadable (the 4c jobs need it).
- A GAP root whose `pkg/primgrp` is this worktree.  It is passed to `-l`, so
  it must keep its trailing `;`.
- Jobs run from the worktree, since the commands name `dev/marta/*.g`.

## Sweep

    marta import-range dev/marta/sweep-range.toml --start 2 --stop 8191

Result is `[degree, groups, complaints]`.  A clean degree has a third entry of
0; anything else means the stored data disagrees with recomputation, and the
complaints are on stdout of the attempt, prefixed `COMPLAINT`.

## Type 4c

    python dev/marta/generate_4c_jobs.py \
        --gaproot '/path/to/gaproot;' \
        --worklist dev/4c-worklist.txt \
        --skip dev/4c-done.txt | marta import -

Result is the wreath decomposition as a flat integer list: `m, k, ngens`, then
per generator the `k` image lists of length `m` and finally sigma as a list of
length `k`.  An entry that cannot be decomposed writes `[0]` and still exits 0
-- a refusal is a real answer here, not a failure, and should not be retried.

Each job verifies its own result before reporting it: the group rebuilt from
the decomposition must equal the original conjugated by the point relabelling.
That is an identity, not an invariant comparison, so a job that reports a
result has proved it.

## Merging results back

The integer lists convert to data-file entries with `PGProductAction4c`
arguments in the obvious way; `dev/product4c.g` has the writer used locally.
Apply them exactly as the local runs do, then re-run the affected degrees of
the sweep, since changing `data/` invalidates any earlier sweep of that degree.

## Known slow spots

- Degrees 4096 and 6561 dominate the sweep.
- Entry 2401/1173 hung locally for hours and is excluded from
  `dev/4c-worklist.txt` callers via `--skip`; give it a longer timeout or
  leave it as permutations.
