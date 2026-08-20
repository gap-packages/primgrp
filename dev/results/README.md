# Computed results not yet applied to data/

These are the outputs of long runs, kept in the repo so a laptop failure or a
handover does not lose hours of computation.  They are inputs to the appliers,
not data files.

| file | what |
|---|---|
| `4c-entries.txt` | `degree nr ["4c",...]` for converted type 4c entries |
| `4c-log.txt` | per-entry outcome and runtime for the 4c run |
| `classical-amb-entries.txt` | the 45 ambiguous type 2 entries resolved by conjugacy-class count |
| `sweep-log.txt` | one line per swept degree; `BAD` lines are complaints |

Applying any of these rewrites `data/`, which invalidates the sweep for every
degree touched.  Delete those degrees' marker files and re-run them.
