#!/bin/sh
#############################################################################
##
##  dev/product4c-run.sh ROOT LOG
##
##  Drive dev/product4c.g until it has converted everything it can.
##
##  A handful of type 4c entries take far longer to split than the rest, and
##  a plain timeout cannot tell one of those from a run that is simply still
##  working.  So this watches the log instead: if nothing new has been printed
##  for STALL seconds, the entry the log last announced is the one to blame.
##  It is killed, added to the skip list, and the run resumed -- the converter
##  rewrites each file as it goes, so nothing but that one entry is lost.
##
##  What ends up in the skip list is the tail worth studying rather than
##  waiting out.
##
ROOT=$1
LOG=$2
STALL=300
ROUNDS=20

skip=""
round=0
while [ $round -lt $ROUNDS ]; do
  round=$((round + 1))
  gap -q -b -A --quitonbreak -l "$ROOT;" \
      -c "conv_dir:=\"data\";; conv_skip:=[$skip];;" dev/product4c.g >> "$LOG" 2>&1 &
  pid=$!

  # kill the run once the log has been quiet for STALL seconds
  while kill -0 $pid 2>/dev/null; do
    sleep 30
    quiet=$(( $(date +%s) - $(stat -f %m "$LOG") ))
    [ $quiet -lt $STALL ] && continue
    kill $pid 2>/dev/null
    wait $pid 2>/dev/null
    stuck=$(grep '^TRY ' "$LOG" | tail -1 | awk '{print "["$2","$3"]"}')
    echo "STALLED $stuck" >> "$LOG"
    skip="$skip${skip:+,}$stuck"
    break
  done

  wait $pid 2>/dev/null && break
done

echo "SKIPPED [$skip]" >> "$LOG"
