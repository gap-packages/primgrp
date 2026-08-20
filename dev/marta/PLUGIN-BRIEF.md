# Brief: a MARTA plugin for primgrp data work

For whoever writes the plugin. This states what the two task families must do,
what they need installed, what they return, and how results come back. It does
not prescribe plugin internals beyond the API v3 surface.

Contact point for questions about the mathematics or the GAP side: the primgrp
branch below carries working standalone versions of both workers, already
tested under MARTA's calling convention. Reuse or rewrite them as you prefer.

## 1. Why

The GAP package `primgrp` ships the primitive groups library. Work on the
branch named below has cut its data from ~2.92 GB to 93.3 MB by replacing
stored permutations with constructions. Two jobs remain, both too slow for one
laptop and both embarrassingly parallel:

- **sweep** — recompute every stored property of every primitive group of
  degree 2..8191 and compare against the stored value. 8190 units.
- **type 4c** — re-encode product-action entries as elements of
  `Sym(m) wr Sym(k)` instead of permutations of `m^k` points. 1553 units,
  currently ~52 MB of the remaining data, expected to become ~1.5 MB.

They are independent of each other and can run concurrently.

## 2. Code and versions

| what | value |
|---|---|
| repository | `git@github.com:gap-packages/primgrp` |
| branch | `claude/data-compression-consistency-8cb77d` |
| commit at time of writing | `d1e0a6a4c3c756474d9c4ae1f41a8328430a54ad` |
| GAP used locally | `4.17dev-75-g525144b` |

Pin an exact commit when you generate jobs and record it. **The data files
change as results come back**, so a result is only meaningful against the
commit that produced it. A sweep result for a degree whose data later changed
must be discarded and the degree re-run; see §8.

GAP packages needed by the two workers:

- `primgrp` itself, from the branch checkout, reachable as a GAP root — that
  is, some directory `R` with `R/pkg/primgrp` being the checkout. Passed to
  GAP as `-l 'R;'`; **the trailing semicolon matters**, it means "and the
  default roots too".
- `json`, for the MARTA GAP shim.
- Nothing else. Both workers were tested with `gap -A` (no autoloading), so
  `transgrp` is *not* required, despite other tooling in `dev/` using it.

`gap` is on `PATH` on tutulla; declare it as a tool anyway so the version is
recorded.

## 3. Task family `primgrp_sweep`

**Input:** one integer, the degree, 2..8191.

**Work:** load `primgrp`, read `tst/testutils.g`, call
`PrimGrpCheckDegree(deg, PrimGrpCheckAll)`. That returns a possibly empty list
of complaint strings, each naming a group and the property that disagreed.

**Output:** the degree, the number of primitive groups of that degree, the
number of complaints, and the complaint strings themselves. A clean degree has
zero complaints — that is the expected case and 5146 degrees have produced it
so far with none failing.

Suggested envelope, following the `grpconst` pattern:

```json
{
  "schema": "marta.primgrp.result.v1",
  "task_family": "primgrp_sweep",
  "degree": 100,
  "result_format": "integer-list",
  "result": [100, 38, 0],
  "complaints": []
}
```

Complaints are free text and may contain quotes and backslashes; if that is
awkward for the result contract, put the count in `result` and the strings in
the attempt's stdout, prefixed `COMPLAINT `. Do not drop them: a complaint is
the entire point of the job.

**Working standalone version:** `dev/marta/sweep-job.g`, driven by
`marta_deg` set via `-c`.

## 4. Task family `primgrp_type4c`

**Input:** two integers, degree and entry number, from
`dev/4c-worklist.txt` (1553 lines, `degree nr`).

**Work:** `Read("dev/product4c.g")`, then `PRIMGRP_Extract4c` on
`PrimitiveGroup(deg, nr)`. This recovers the identification of the `m^k`
points with `k`-tuples and rewrites each generator as `(p_1,...,p_k; sigma)`.

**Output:** `m`, `k`, and per generator the `k` permutations of degree `m`
plus `sigma` of degree `k`, each as an integer image list.

```json
{
  "schema": "marta.primgrp.result.v1",
  "task_family": "primgrp_type4c",
  "degree": 25, "nr": 24,
  "result_format": "product-action-decomposition/v1",
  "result": {
    "m": 5, "k": 2,
    "gens": [
      {"base": [[1,2,4,5,3],[2,3,4,5,1]], "top": [1,2]}
    ]
  }
}
```

`permutation-group-generators/v1` does not fit: these are wreath elements, not
generators of a permutation group on the job's own degree. A new format is
warranted. A flat integer list `[m, k, ngens, base..., top...]` also works if
you prefer to reuse `integer-list`; `dev/marta/4c-job.g` currently emits that.

**Self-verification is mandatory and already implemented.** Before reporting,
the worker rebuilds the group from the decomposition and checks it equals the
original conjugated by the point relabelling. That is an identity, not a
comparison of invariants, so a reported result is proved correct and the merge
step needs no further checking. Do not remove this. It is cheap relative to
the extraction.

**Working standalone version:** `dev/marta/4c-job.g`, driven by `marta_deg`
and `marta_nr`.

## 5. Refusal is not failure

Some units cannot be done, for mathematical reasons, and will never succeed on
retry:

- 4c entries whose socle does not decompose as expected.
- One known unit, degree 2401 entry 1173, hung locally for 5.8 hours and is
  excluded via `--skip`. Give it a long timeout in a one-off run if you want an
  answer, but do not let it block a batch.

Both workers signal this by returning a refusal **with exit status 0** and a
reason on stdout prefixed `REFUSED `. `max_attempts` must not cause a retry of
a refusal — retrying only spends the time again. If the plugin distinguishes
"succeeded with no result" from "failed", use that; otherwise a result of
`[0]` / `{"m": 0}` means refused.

Genuine failures — GAP crash, out of memory, timeout — should retry.

## 6. Resource profile, measured

Per-unit wall clock on one laptop core:

| | n | median | p90 | p99 | max |
|---|---|---|---|---|---|
| sweep | 5218 | 9 ms | 39 s | 209 s | 2.55 h (degree 4096) |
| 4c | 809 | 293 ms | 21 s | 679 s | 66 min |

Two consequences:

- **One job per unit, never per block.** The distribution is extreme — half
  the sweep degrees finish in under 10 ms while degree 4096 takes 2.55 hours.
  Blocking leaves one worker running for hours while the rest idle, which is
  what happened locally.
- **Timeouts must be generous.** A 10-minute cap would lose the top 1% of both
  families, which is precisely the expensive, interesting part. Suggest 4 h for
  sweep and 2 h for 4c, or no cap with a slot budget.

Memory: a long-lived GAP accumulated to 450 MB and was killed by the OS twice
locally, because `PRIMGrp` memoises every entry it evaluates. **One GAP process
per unit avoids this entirely** — per-unit footprint is small. If you batch
units into one process for startup amortisation, cap the batch and restart.

GAP startup with `-A` is ~1 s; against a 293 ms median for 4c that is real
overhead. Batching 4c units by degree would amortise it and share the
degree's loaded data file, at the cost of reintroducing the memory issue.
Batching within a degree is the safer compromise if you want it.

## 7. Progress already made

Do not redo work. Skip lists are in the branch:

| file | lines | meaning |
|---|---|---|
| `dev/sweep-done.txt` | 5146 | degrees already swept clean at commit `d1e0a6a` |
| `dev/4c-done.txt` | 809 | 4c entries already converted |
| `dev/4c-worklist.txt` | 1553 | all 4c units |

`dev/marta/generate_4c_jobs.py` takes `--worklist` and `--skip` and emits
JSONL for `marta import -`. It also takes `--count`.

**Caveat on `dev/sweep-done.txt`:** those degrees were swept against data that
has since been partly rewritten. Any degree whose entries changed after it was
swept must be re-run. If in doubt, re-run the whole sweep; it is cheap in
aggregate once parallel.

## 8. Getting results back

Results must end up as edits to `data/gps*.g` on the branch. Proposed flow,
but do whatever suits MARTA:

1. Jobs write results into MARTA as above.
2. An export command emits `degree nr <entry-text>` lines, the format the
   local appliers already consume.
3. Those are applied to `data/`, the test suite is run, and the change is
   committed.
4. Every degree touched by an applied change is re-run through the sweep,
   because changing `data/` invalidates any earlier sweep of that degree.

Step 4 is the one that is easy to forget and the reason results are tied to a
commit. Locally this is handled by deleting the marker files of affected
degrees.

Transfer: `git clone` of the branch to tutulla is enough for the code; the
data files are in the repo, so nothing extra to rsync in that direction.
Coming back, the entry texts are small — the whole 4c output is expected to be
about 1.5 MB — so they can travel through MARTA's artifact store rather than
by `scp`.

## 9. What I will check when reviewing

- The GAP root is passed with its trailing `;` and jobs run from the checkout.
- 4c self-verification is intact and its result is not trusted without it.
- Refusals do not retry; genuine failures do.
- Timeouts do not truncate the expensive tail.
- One GAP process per unit, or a bounded batch with restarts.
- The commit that produced each result is recorded, and the sweep re-run rule
  in §8 is expressible.
- A smoke run of one unit from each family reproduces the local results:
  degree 100 gives `[100, 38, 0]`, and degree 25 entry 24 gives `m=5, k=2`
  with 4 generators.
