# Running the long primgrp jobs elsewhere

Two jobs dominate this branch and neither fits comfortably on a laptop: the
consistency sweep over degrees 2..8191, and the type 4c re-encoding of 1553
entries.  Both are embarrassingly parallel and resumable, so they shard.

## What has to be true on the target

- `gap` on PATH, with `transgrp` available (the 4c work needs it).
- This worktree reachable as a GAP root: `PRIMGRP_ROOT=/some/dir` where
  `$PRIMGRP_ROOT/pkg/primgrp` is (or symlinks to) the worktree.
- A shared or collected output directory.

## Submitting

Slurm, 32 shards:

    sbatch --array=0-31 --wrap \
      "dev/cluster/shard.sh sweep \$SLURM_ARRAY_TASK_ID 32 \$SCRATCH/sweep"

    sbatch --array=0-31 --wrap \
      "dev/cluster/shard.sh 4c \$SLURM_ARRAY_TASK_ID 32 \$SCRATCH/4c dev/4c-worklist.txt"

Without a scheduler, the same script takes the shard number as an argument, so
`for i in $(seq 0 31); do dev/cluster/shard.sh 4c $i 32 out & done` works.

## Cost

Measured on one laptop core: the sweep reached 4869 of 8190 degrees in 20.5 h,
with degree 4096 alone taking 2.55 h; 4c did 729 of 1553 entries in 13 h, but
7.2 h of that was one entry that hung.  Sharding 32 ways should bring both
under an hour of wall clock, with the caveat that the cost per degree is very
uneven -- hence the striping in the sweep shards rather than contiguous blocks.

## Collecting

Each shard writes `entries-<shard>.txt` and `log-<shard>.txt`.  Concatenate
the entries files and apply them exactly as the local runs do; the marker
directory can be shared or merged afterwards.

## Known slow spots

- degree 4096 and 6561 dominate the sweep.
- 2401 nr 1173 hung for hours locally and is worth a longer TIMEOUT or a skip.
