# Coarse forward guide: numerical screen passed, sampling unproven

Preflight job `122252` completed successfully. It used three patients (3,74,122),
all eight local coordinates, three dynamic-shared contexts and three fixed
tolerance profiles. The exact target, observation model and patient time grids
were unchanged. Thirty implementation assertions passed before the preflight.

| VODE-BDF guide relative tolerance | Median CPU speedup | Smallest case speedup | Largest case 95% absolute residual difference | Preflight go |
|---|---:|---:|---:|---|
| 0.0001 | 1.574 | 1.538 | 0.00241 | No: insufficient median speed |
| 0.01 | 2.596 | 2.174 | 0.0935 | Yes |
| 0.1 | 3.099 | 2.691 | 0.935 | Yes |

All exact design evaluations were finite, and there were no approximate solver
failures on those states. Fresh repeated likelihood evaluations were identical.
There were 1,728 exact and 5,184 guide calls; diagnostic elapsed time was 34.6
seconds. Complete Slurm preflight time was 46 seconds and TotalCPU was 40.182
seconds, including setup/tests. That complete CPU cost is charged to the guide
in the following comparison, including profiles that were not selected.

Following the [frozen preflight rule](COARSE_GUIDE_CONTRACT.md), select the
fastest qualifying profile, relative tolerance 0.1 and absolute tolerance 0.001.
The residual figures compare disjoint pairs of proposal-design states. They
are NOT stationary acceptance or mixing estimates. A 3.1-fold cheaper numerical
guide is not a 3.1-fold faster sampler.

## Actual sampling test

The [conditional-chain contract](COARSE_GUIDE_CHAIN_CONTRACT.md) compares a
preconditioned exact random walk with eight cheap inner moves followed by an
exact endpoint correction. Four seeds per method, patient and fixed context
give 72 short chains on three patients. The baseline receives 6,000 transitions;
the macro arm 1,500. Actual costs and population-score mixing decide the result.
Implementation smoke `122253` passed the full test suite and short real-target
baseline/macro runs; its 33 exact and 65 guide calls are implementation overhead,
not scientific performance evidence. Full comparison job `122254` is running
on four CPU workers. Its pre-run full test suite passed and initial conditional
chains are completing without errors. No retained comparison is available yet.

This is only an enabling experiment for modular Bayes. A pass would need to be
followed by posterior-accurate reuse across changed population structures with
shared uncertainty restored and all incremental ODE costs counted. It does not
authorise a 115-patient run, a message/evidence claim, or a Gaussian restriction
on future scientific population distributions. No System B source is modified.

[Compact preflight evidence](results/coarse_guide) is committed; full call
records remain in `outputs/coarse_guide_122252/`.
