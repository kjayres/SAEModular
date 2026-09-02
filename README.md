# SAEModular

This research prototype asks a narrow question: can patient transports learned
around a SAEM-like anchor make exact global Bayesian MCMC substantially faster
without changing the posterior target?

For frozen patient maps, the implemented proposal is

```text
u_i     = L_i(c)^(-1) (x_i - m_i(c))
x_i_new = m_i(c_new) + L_i(c_new) u_i
```

It is corrected against the original joint density with the exact affine
Jacobian. Approximate maps therefore affect efficiency, not the invariant
posterior. The move must be combined with an exact local-patient refresh because
it preserves every `u_i` and is not irreducible alone.

## Current evidence

The nonlinear decay toy passed a frozen robustness run against a high-accuracy
quadrature posterior. Affine transport gave:

- 2.676 times the fixed-state ESS per patient-likelihood call;
- 1.780 times the fixed-state posterior-scaled ESJD;
- split R-hat 1.0004;
- posterior-mean error 0.0084 reference standard deviations.

This validates the exact transport mechanism on the toy only. The affine tuning
optimum was at the largest tested proposal scale, so the efficiency comparison
is not claimed to be fully optimized. ESS is the script's transparent
initial-positive-sequence estimate and is corroborated by ESJD and quadrature
agreement.

The System A experiment is intentionally limited to the 12 previously audited
patients. It compares fixed patient states, conditional-mean translation, and
full eight-dimensional affine transport for identical dynamic-shared-parameter
endpoints. It reports
exact VODE-BDF acceptance, whitened ESJD per proposed solve, solver failures,
patient and aggregate work variance, and forward/reverse D2 diagnostics.

The first short-bank pilot did **not** pass: full affine transport achieved only
0.538 times the fixed-state ESJD per ODE at a one-standard-deviation endpoint.
Its conditional banks also failed their independent quality gate (minimum
coordinate ESS 4.4), and the original interpolation did not reproduce its own
training endpoints.

The bounded rescue corrected those defects and produced adequate banks. It
still failed the predeclared 1.5-fold continuation gate: mean-only transport
gave 1.366 times fixed-state ESJD per ODE and full affine transport gave 1.107.
The mean effect was concentrated in the second whitened shared-parameter
direction; covariance transport added no value. Full-cohort and retained-chain
work is therefore stopped. See `RESULTS.md`.

That System A experiment is an **oracle-anchor endpoint pilot**. It is not yet a
SAEM handoff, retained posterior chain, bridge-sampling method, population
parameter sampler, full-115 scaling result, or deployable proof. Its acceptance
calculation assumes the audited symmetric additive proposal in dynamic `psi`.

## Layout

```text
R/affine_map.R                 frozen affine-map primitives
R/transport_mh.R               exact deterministic-transport MH kernel
R/system_a_adapter.R           fail-closed adapter to the certified target
experiments/toy_decay.R        nonlinear toy validation
experiments/system_a_12.R      12-patient oracle endpoint pilot
slurm/*.sbatch                 Oxford Statistics HPC launchers
tests/testthat/                inverse, Jacobian, cache, and target checks
```

All numerical experiments run through Slurm; no ODE work is run on the login
node.

```bash
sbatch slurm/toy_decay.sbatch
sbatch slurm/system_a_12.sbatch
```

Generated logs and RDS/CSV outputs are deliberately excluded from Git.

## System A prerequisites

The System A pilot is not standalone. It requires the local `modular_bayes`
target and comparator artifacts, including the clean pinned commit
`a3e0367b06c82ad4d07280c99f10d4f9bac69978`. The adapter checks commit identity,
source hashes, target fingerprint, and the numerical-audit certificate before
loading the sealed VODE-BDF callbacks. The Slurm launchers contain
Oxford/user-specific absolute paths and must be adapted elsewhere.
