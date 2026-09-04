# SAEModular: fast hierarchical inference experiments

The objective is useful Bayesian population uncertainty and reusable expensive
forward calculations. Engineering is subordinate to demonstrating computational
and statistical benefit. Nothing here establishes a general modular Bayes solution.

## Current bounded experiment

Can a continuously refreshed pool reuse trajectories across System A patients?
The [frozen feasibility contract](SHARED_POOL_CONTRACT.md) asks only whether
candidates are useful to multiple patients and whether a common-trajectory
interface is accurate and cheaper. No full sampler, SAEM rerun, learned density
or neural surrogate is implemented. The sealed target remains unchanged.

**Result:** promising trajectory reuse, but a narrowly failed numerical gate;
no automatic progression to retained chains. Read the
[result and interpretation](SHARED_POOL_RESULT.md) and
[compact diagnostic CSVs](results/shared_pool_probe).

Code is deliberately limited to:

- `R/system_a_adapter.R`: pinned original target and provenance checks;
- `R/shared_pool_snapshots.R`: two coherent Bayesian comparator snapshots;
- `R/shared_trajectory_probe.R`: experimental union-grid versus original solves;
- `R/shared_pool_probe.R`: exact proposal probabilities and offline summaries;
- `experiments/system_a_shared_pool_probe.R`: one bounded experiment driver.

All tests and numerical work run through Slurm, never on the login node:

```sh
sbatch --parsable slurm/test_core.sbatch
sbatch --parsable slurm/system_a_shared_pool_probe.sbatch
```

The experiment driver also accepts `--replay-saved <completed-source-dir>` before
its new output directory. This recomputes summaries from immutable saved endpoint
evaluations without any new scientific ODE calls; it verifies source hashes and
identical contexts/reserves. It must also be run within a Slurm job.

## Archived Fable work

The affine-map, defensive-mixture and pure-SAEM-message implementations have
been removed from this branch. Complete code is recoverable at Git tag
`archive/fable-final-20260904` (commit `bb50aee`). Historical findings remain in
[RESULTS.md](RESULTS.md), [REVIEW_PACKET.md](REVIEW_PACKET.md), and the
[final pure-SAEM result](FINAL_PURE_SAEM_RESULT.md) and
[contract](FINAL_PURE_SAEM_CONTRACT.md). References to old source paths in those
reports refer to that archived revision.

This retirement is a computational decision about the tested implementations,
not a proof that all SAEM-conditioned messages or all patient samplers fail.
The likelihood-based SAEM endpoints also did not incorporate the informative
population hyperpriors of the Bayesian target.

Raw research outputs remain ignored by Git. Compact evidence belongs in the
repository; expensive-to-regenerate banks are not deleted simply to save a few MB.
