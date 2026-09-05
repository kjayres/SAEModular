# SAEModular: fast hierarchical inference experiments

The objective is useful Bayesian population uncertainty and reusable expensive
forward calculations. Engineering is subordinate to demonstrating computational
and statistical benefit. Nothing here establishes a general modular Bayes solution.

## Current reduced posterior experiment

An exact-corrected shared-pool sampler is being tested against the same baseline
hierarchy on twelve patients, with four chains per method and a fixed total
evaluation budget. The original sealed likelihood and informative priors define
both targets; common trajectories only guide proposals. There are no retained
results yet. The [frozen chain contract](SHARED_POOL_CHAIN_CONTRACT.md) requires
posterior agreement and a material improvement in population ESS per total CPU
cost before further work. It includes dynamic shared-parameter uncertainty.

## Residual-collapse prerequisite: checked, sampler parked

The [one-patient bound audit](RESIDUAL_BOUND_RESULT.md) used 26 sealed prediction
calls. The proposed global Gaussian/affine lower-bound construction does not
apply unchanged to System A. Simple Gaussian-integrable regions repair physical
equilibrium support only; no useful certified likelihood lower component has
been constructed. The general residual identity is not falsified. No residual
sampler or certification framework was built. See the
[frozen scope](RESIDUAL_BOUND_FEASIBILITY.md) and
[compact evidence](results/residual_bound_probe).

This repository is the bounded methods testbed. Promote validated components to
the main modular Bayes work only after the relevant posterior-accuracy and total-
cost tests; a sampling optimisation is not automatically a modularisation result.
Retain compact negative evidence and distinguish failed constructions from
inconclusive tests or general impossibility claims.

## Completed feasibility screen

Can a continuously refreshed pool reuse trajectories across System A patients?
The [frozen feasibility contract](SHARED_POOL_CONTRACT.md) asks only whether
candidates are useful to multiple patients and whether a common-trajectory
interface is accurate and cheaper. That screen did not run a retained sampler,
SAEM rerun, learned density or neural surrogate. The sealed target remains unchanged.

**Result:** promising trajectory reuse, but a narrowly failed numerical gate;
no automatic progression to retained chains. Read the
[result and interpretation](SHARED_POOL_RESULT.md) and
[compact diagnostic CSVs](results/shared_pool_probe).

Main code:

- `R/system_a_adapter.R`: pinned original target and provenance checks;
- `R/shared_pool_snapshots.R`: hash-pinned coherent Bayesian comparator snapshots;
- `R/shared_trajectory_probe.R`: experimental union-grid versus original solves;
- `R/shared_pool_probe.R`: exact proposal probabilities and offline summaries;
- `experiments/system_a_shared_pool_probe.R`: one bounded experiment driver.
- `R/shared_pool_kernel.R`, `R/system_a_population_updates.R`: corrected
  assignment and exact population kernels;
- `R/system_a_shared_chain.R`: the reduced baseline/recycling chain;
- `experiments/system_a_shared_pool_chain.R`: paired runs, cost and posterior diagnostics.

All tests and numerical work run through Slurm, never on the login node:

```sh
sbatch --parsable slurm/test_core.sbatch
sbatch --parsable slurm/system_a_shared_pool_probe.sbatch
sbatch --parsable --export=ALL,SMOKE=1 slurm/system_a_shared_pool_chain.sbatch
sbatch --parsable slurm/system_a_shared_pool_chain.sbatch
sbatch --parsable slurm/system_a_residual_bound_probe.sbatch
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
