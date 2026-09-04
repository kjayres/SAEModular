# Shared-trajectory feasibility result

2026-09-05 (Europe/London). Scientific job 122228, source `4963c95`;
corrected, ODE-free reanalysis job 122229, source `c4746e8`. Both jobs completed
with exit 0. The full original suite passed; the trimmed suite and new timing
regression also passed. No SAEM refit, conditional-bank simulation, retained
sampler, neural training or population-structure comparison was run.

## Decision

**Evidence of a useful computational ingredient, not a demonstrated fast
Bayesian sampler. The predeclared overall screen did not pass.**

There is genuine trajectory reuse, and the common-grid calculation is materially
cheaper. The sole frozen gate failure is numerical disagreement at two extremely
implausible proposed patient states. Keep that failure recorded; do not change
the threshold or call the new numerical target certified. Equally, this is not
evidence that the statistical recycling mechanism failed.

## The two questions

Both contexts are coherent joint Bayesian comparator draws (chains 1 and 3,
retained draw 750), not the likelihood-only SAEM endpoints. The same 128 direct
draws from a fixed four-population mixture are evaluated at both contexts.

| Measurement | Context 1 | Context 2 |
|---|---:|---:|
| Common evaluation speedup, 12-patient panel | 7.16x | 7.03x |
| Common evaluation speedup, 115-patient grid audit | 26.61x | 22.53x |
| Reserves useful to >=2 panel patients | 16/128 (12.5%) | 46/128 (35.9%) |
| Reserves useful to >=2 **ODE-requiring** patients | 8/128 (6.25%) | 17/128 (13.28%) |
| Mean one-reserve acceptance, all panel patients | 2.26% | 5.23% |
| Patients with mean one-reserve acceptance >=2% | 6/12 | 7/12 |
| Largest absolute log-likelihood difference | 0.012811 | 0.007160 |
| Largest change in reserve MH probability | 0.00000991 | 0.00001943 |
| Finite-support / prediction-outcome disagreements | 0 / 0 | 0 / 0 |
| Frozen numerical gate (maximum error <=0.01) | fail | pass |

“Useful” means one-reserve MH probability >=0.1. These are potential alternative
uses of a candidate, not multiple simultaneous occupants or achieved chain moves.
The timing ratio compares evaluation of the SAME vector for every patient with
one union-grid solve plus patient likelihoods. It is NOT a sampler speedup and
does not include a future live pool's assignment bookkeeping/cache invalidation.
The full-grid audit contains only eight predeclared vectors per context.

Patients 20, 71 and 117 need equilibrium calculations only. The original
all-patient usefulness gate counted them; the ODE-only breakdown was added during
review and does not replace the frozen gate. It is the more relevant description
of expensive cross-patient reuse. Some apparent sharing was already ODE-free.

Coverage is uneven. Patient 6 has direct-reserve acceptance of only 0.0276% and
0.0788%, and no reserve crosses the 0.1 usefulness threshold for that patient.
Patients 3 and 55 also have weak direct-reserve coverage. The pool cannot yet be
treated as a replacement for all patient-local updates. Existing patient states
have mean cross-patient acceptance 5.11% and 9.92% excluding self-comparisons,
but those states must first be vacated before they can legally be reserves.

## Why the numerical failure does not condemn recycling

The two threshold violations are patient 3 at reserves 28 and 46 in context 1.
Original log likelihoods are -2796.01 and -2985.38; the current state's is about
-59.18. Both exact and common-grid proposal acceptances underflow to zero.
The [two raw rows](results/shared_pool_probe/numerical_gate_failures.csv) are
preserved. Their absolute errors breach the declared limit, but they do not
meaningfully change these particular proposal decisions.

The interface remains an experimental numerical target. Agreement on these
states is not proof of agreement or matching solver-failure boundaries everywhere.

## Costs and reporting correction

- Scientific prediction-interface calls: 5,496 (exact predeclared cap).
- Actual integration attempts across both interfaces: 4,566.
- Integration failures: zero. Pre-integration equilibrium rejection remains
  counted in candidate denominators; equilibrium-only evaluations are not ODEs.
- Scientific driver wall time, including loading and reporting: 92.692 seconds
  on four allocated CPUs. Entire job including the old regression suite: 9m21s.
- Reanalysis 122229: **zero new scientific ODE calls**, entire job 2m54s.

The original timing CSV incorrectly pooled both grid sizes because R's
`subset()` masked the loop variable with its `grid` column. The raw per-vector
ledger was correct. Explicit grouping and a regression test repaired the report
from saved, hash-verified evaluations. The corrected per-grid speed gates both
still pass. No thresholds, proposals, contexts or numerical evaluations changed.
Original outputs stay in `outputs/shared_pool_probe_122228`; corrected reports
are in `outputs/shared_pool_probe_122229`.

Compact, corrected evidence is committed in
[results/shared_pool_probe](results/shared_pool_probe): decisions, timings,
patient and reserve summaries, complete call ledger, source hashes and the two
numerical gate failures. Raw states and observations are not committed.

## What this suggests next, if separately authorized

The present uncorrected common-grid route stops at its failed numerical gate;
no retained chain follows automatically. There is nevertheless a concrete way to
use the promising part without silently changing the target: use common-grid
outputs only to select cheap proposals, then apply an exact correction using
the original patient likelihood. This is the same reversible-surrogate principle
as the proposed learned-forward-model approach, but requires no neural training
to obtain a first cheap approximation. A poor approximation affects efficiency,
not the corrected invariant target, when the proposal/support conditions hold.
See [Multilevel Delayed Acceptance MCMC](https://epubs.siam.org/doi/10.1137/22M1476770).

That is a proposed new kernel, not something this screen has implemented or
validated. Its next scientific test would be global posterior accuracy and ESS
per total cost in a reduced hierarchy, including reserve turnover, exact
corrections, ordinary local moves and dynamic-shared-parameter cache rebuilding.
These measurements, rather than more auxiliary fitting, must decide whether it
is worth pursuing. Neither this screen nor the older failures establish a
general solution or impossibility theorem for modular Bayes.

## Repository cleanup

About 14,000 lines of retired transport/SAEM-message code, tests and launchers
were removed from this branch; all remain recoverable at
`archive/fable-final-20260904` (`bb50aee`). Historical scientific reports and
expensive-to-regenerate raw results were retained. Two verbose old SAEM failure
JSON logs were losslessly compressed from 12,760,001 to 1,032,708 bytes, freeing
about 11.2 MiB; `gunzip` restores their original contents and filenames. One
untracked default `Rplots.pdf` scratch plot was discarded (not Git-recoverable).
No external project, source dataset or pinned target worktree was changed.
