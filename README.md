# SAEModular

This research prototype tests one idea: use an all-115-patient SAEM fit as a
hierarchy-aware anchor, then turn independent exact patient-conditional banks
into reusable messages for fast modular population inference. At fixed shared
parameters, patient `i` contributes the raw-bank message

```text
r_i(M, eta) = mean_s g_(i,M)(x_is | eta) / g_(i,A)(x_is | eta_A),
x_is ~ p(x_i | y_i, eta_A, psi_A).
```

Once the exact conditional banks have been built, evaluating these messages is
ODE-free. Raw-bank importance reweighting is tested first. Quadratic surfaces,
transport maps, changing dynamic shared parameters, and retained full-cohort
chains are deliberately outside the first test.

The first falsification sequence is:

1. fit the genuine all-patient SAEM anchor;
2. freeze all four shared parameters at that anchor;
3. build independent exact VODE-BDF conditional banks for the 12 audited
   patients;
4. calculate raw `g_new / g_anchor` messages without ODE calls;
5. validate selected population endpoints with independent candidate banks and
   bidirectional bridge estimates;
6. proceed to all 115 patients only if overlap, replication, and accuracy pass.

The frozen-bank posterior is a controlled Monte Carlo approximation, not an
exact retained MCMC target. Independent banks, weight diagnostics, and bridge
checks quantify whether that approximation is accurate enough. Dynamic shared
parameter uncertainty is a later problem because changing those parameters
requires new likelihood/ODE evaluations.

## Earlier transport result

Before the message-first scope was clarified, this repository tested an exact
affine global co-move. The nonlinear decay toy passed, but the corrected
12-patient System A rescue missed its predeclared continuation gate: mean-only
transport achieved 1.366 times fixed-state ESJD per ODE and full affine
transport 1.107, below the required 1.5. That branch is stopped and retained
only as negative evidence; see `RESULTS.md`.

## Layout

```text
R/affine_map.R                 frozen affine-map primitives
R/transport_mh.R               exact deterministic-transport MH kernel
R/system_a_adapter.R           fail-closed adapter to the certified target
R/system_a_patient_bank.R      exact fixed-anchor patient pCN banks
R/patient_messages.R           raw-message and bridge estimators
experiments/toy_decay.R        nonlinear toy validation
experiments/system_a_12.R      12-patient oracle endpoint pilot
config/system_a_saem_sources.sha256  pinned external SAEM inputs
slurm/system_a_saem_anchor.sbatch    all-115 anchor-only SAEM fit
slurm/test_core.sbatch              compute-node source/synthetic tests
slurm/*.sbatch                 Oxford Statistics HPC launchers
tests/testthat/                inverse, Jacobian, cache, and target checks
```

All numerical experiments run through Slurm; no ODE work is run on the login
node.

```bash
sbatch slurm/toy_decay.sbatch
sbatch slurm/system_a_12.sbatch
sbatch slurm/system_a_saem_anchor.sbatch
sbatch slurm/test_core.sbatch
```

Generated logs and RDS/CSV outputs are deliberately excluded from Git.

## System A prerequisites

The System A pilot is not standalone. It requires the local `modular_bayes`
target and comparator artifacts, including the clean pinned commit
`a3e0367b06c82ad4d07280c99f10d4f9bac69978`. The adapter checks commit identity,
source hashes, target fingerprint, and the numerical-audit certificate before
loading the sealed VODE-BDF callbacks. The Slurm launchers contain
Oxford/user-specific absolute paths and must be adapted elsewhere.
