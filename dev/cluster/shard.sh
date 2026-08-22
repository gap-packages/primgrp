#!/bin/bash
# Run one shard of a long primgrp job.  Both jobs below are embarrassingly
# parallel and resumable: every unit writes a marker file when attempted, so
# shards never collide and a re-submission picks up where it stopped.
#
#   ./shard.sh sweep <shard> <nshards> <outdir>
#   ./shard.sh 4c    <shard> <nshards> <outdir> <worklist>
#
# <shard> comes from the scheduler if not given: SLURM_ARRAY_TASK_ID,
# SGE_TASK_ID or PBS_ARRAYID.  GAP must be on PATH and the primgrp worktree
# reachable as a GAP root (see ROOT below).
set -u
job=${1:?job: sweep|4c}
shard=${2:-${SLURM_ARRAY_TASK_ID:-${SGE_TASK_ID:-${PBS_ARRAYID:-0}}}}
n=${3:?number of shards}
out=${4:?output directory}
ROOT=${PRIMGRP_ROOT:?set PRIMGRP_ROOT to a dir whose pkg/primgrp is the worktree}
CHUNK=${CHUNK:-60}

mkdir -p "$out/marks"
cd "$(dirname "$0")/../.."

case "$job" in
sweep)
  # degrees 2..8191, striped so every shard gets a mix of cheap and dear ones
  for ((lo=2+shard; lo<=8191; lo+=n*CHUNK)); do
    hi=$((lo + (n*CHUNK) - 1)); ((hi > 8191)) && hi=8191
    cat > "$out/chunk-$shard.g" <<GEOF
PRIMGRP_Degrees := Filtered([$lo..$hi], d -> d mod $n = $shard mod $n);;
PRIMGRP_MarkDir := "$out/marks";;
PRIMGRP_LogFile := "$out/log-$shard.txt";;
Read("dev/check-sweep.g");
GEOF
    timeout "${TIMEOUT:-7200}" gap -q -b -A --quitonbreak -l "$ROOT;" "$out/chunk-$shard.g" \
        >> "$out/run-$shard.log" 2>&1
  done
  ;;
4c)
  list=${5:?worklist file}
  awk -v s="$shard" -v n="$n" 'NR % n == s % n' "$list" > "$out/work-$shard.txt"
  split -l "$CHUNK" "$out/work-$shard.txt" "$out/part-$shard-"
  for part in "$out/part-$shard-"*; do
    cat > "$out/chunk-$shard.g" <<GEOF
Read("dev/product4c.g");
BreakOnError := false;;
out := OutputTextFile("$out/entries-$shard.txt", true);; SetPrintFormattingStatus(out,false);;
log := OutputTextFile("$out/log-$shard.txt", true);; SetPrintFormattingStatus(log,false);;
for l in SplitString(StringFile("$part"), "\n") do
  if Length(l) = 0 then continue; fi;
  p := List(SplitString(l, " "), Int);
  mark := Concatenation("$out/marks/", String(p[1]), "_", String(p[2]));
  if IsReadableFile(mark) then continue; fi;
  FileString(mark, "");          # before attempting: a hang costs one chunk
  t := Runtime();
  r := CALL_WITH_CATCH(PRIMGRP_Convert4c, [p[1], p[2]]);
  if not r[1] then PrintTo(log, p[1], " ", p[2], " EXCEPTION\n");
  elif IsBound(r[2].err) then PrintTo(log, p[1], " ", p[2], " SKIP ", r[2].err, "\n");
  else
    PrintTo(out, p[1], " ", p[2], " ", r[2].text, "\n");
    PrintTo(log, p[1], " ", p[2], " OK ", Length(r[2].text), " ms ", Runtime()-t, "\n");
  fi;
od;
CloseStream(out); CloseStream(log);
QUIT;
GEOF
    timeout "${TIMEOUT:-7200}" gap -q -b -A --quitonbreak -l "$ROOT;" "$out/chunk-$shard.g" \
        >> "$out/run-$shard.log" 2>&1
  done
  ;;
*) echo "unknown job $job" >&2; exit 1;;
esac
