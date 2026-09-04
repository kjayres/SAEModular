# SAEModular

This research prototype tests one idea: use an all-115-patient SAEM fit as a
hierarchy-aware anchor, then turn independent exact-target patient-conditional
MCMC banks into reusable messages for fast modular population inference. The
current and final System A pilot samples directly under the fitted SAEM
population distribution:

```text
h_i(x) = g_SAEM,i(x).
```

With the patient likelihood unchanged and all shared parameters fixed, patient
`i` contributes the raw-bank message

```text
r_i(M, eta) = mean_s g_(i,M)(x_is | eta, z_i) / h_i(x_is),
x_is ~ L_i(y_i | x_i, psi_A) h_i(x_i).
```

Once the conditional MCMC banks have been built, evaluating these messages is
ODE-free. The banks target the exact likelihood but finite warm-up does not make
their retained states iid exact conditional draws; convergence and replication
remain experimental requirements. Raw-bank importance reweighting is tested
first. Quadratic surfaces, transport maps, changing dynamic shared parameters,
and retained full-cohort chains are deliberately outside the first test.

The message identity makes no Gaussian or diagonal-covariance assumption.
`g_(i,M)` may be a correlated distribution, heavy-tailed family, mixture, or a
population model with patient covariates `z_i` and discrete structure `M`, as
long as its normalized log density can be evaluated on every stored local
state. Conditional on `(M, eta)`, the population law must factor into these
patient-specific densities; cross-patient latent dependence would require a
larger joint module rather than a product of patient messages. Model-dependent
normalizing constants must be included. Candidate and anchor models must use
the same patient likelihood, local-state coordinates, and dominating measure,
and the anchor support must cover the candidate's likelihood-weighted support.
A covariate may enter `g_(i,M)` without an ODE call. If a changed covariate or
model also changes the observation likelihood or ODE dynamics, however, that
likelihood ratio does not cancel and the move is not ODE-free.
A model that changes the likelihood, dimension, or interpretation of `x_i`
needs an explicit augmentation/mapping, likelihood ratio, or separate bank;
the population-only ratio above is then not valid by itself.

Support and tail adequacy are caller obligations and cannot be established from
a finite bank. The candidate likelihood-weighted measure must be absolutely
continuous with respect to the anchor measure and have a finite marginal
likelihood. Useful Monte Carlo accuracy additionally requires adequate overlap;
standard MCSE arguments require suitable weight moments. In particular, a
Gaussian anchor can have full support yet give infinite-variance or
catastrophically rare weights for a heavier-tailed candidate. Empirical ESS,
split checks, and D2 can expose observed failure but cannot rule out unseen tail
failure.

The final System A implementation is narrower than this generic identity. Its
pure SAEM reference is an eight-dimensional diagonal Gaussian whose mean
includes the observed Nelfinavir indicator. A pCN proposal preserves that
reference, so its exact MH correction is only the sealed likelihood ratio. It
does **not** require later candidate message models `g_(i,M)` to be diagonal
Gaussian: those models enter through normalized log-density evaluations. An
arbitrary non-Gaussian anchor would require a different exact-target
conditional sampler, and bidirectional validation requires a valid sampler at
the candidate target.

The earlier equal mixture
`0.5*g_SAEM + 0.5*g_priorcentral` was a deliberately broad defensive anchor.
It mixed poorly and failed its predeclared calibration, so that branch is
stopped; it remains in the repository only as negative evidence.

The repository now exposes the distribution-free log-density-ratio estimator,
but it does not yet contain a catalogue of candidate `g_(i,M)` models or the
outer sampler over `(M, eta)`. The first message experiment must supply those
normalized candidate log densities explicitly and keep their covariate data
and model priors under provenance control.

The final bounded falsification sequence is:

1. use both corrected all-patient SAEM endpoints and freeze each branch's four
   shared parameters;
2. build four independent pure-SAEM exact-target pCN chains for each of the 12
   audited patients under a fixed short ODE budget;
3. use chains 01-02 alone to test four meaningful population changes: a
   one-population-SD treatment shift in both directions and a 0.8/1.25 change
   in the lambda population SD;
4. stop before candidate ODE work unless every endpoint passes at both SAEM
   branches;
5. if both pass, use held-out reference chains, separately simulated candidate
   banks, and disjoint raw/bridge chain roles to compare forward, reverse, and
   bridge messages;
6. license a 115-patient fixed-shared-parameter test only if all final gates
   pass. The frozen contract is in `FINAL_PURE_SAEM_CONTRACT.md`.

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
R/system_a_patient_bank.R      System A-specific diagonal-prior pCN banks
R/patient_messages.R           distribution-free raw-message/bridge estimators
R/system_a_saem_anchor.R       strict SAEM-to-System-A anchor conversion
R/system_a_message_validation.R defensive-reference validation
R/system_a_pure_message_validation.R final pure-SAEM validation and gates
experiments/toy_decay.R        nonlinear toy validation
experiments/system_a_12.R      12-patient oracle endpoint pilot
experiments/system_a_anchor_bank_patient.R  one-patient exact-target MCMC worker
experiments/system_a_pure_message_validation.R final plan/assessment CLI
config/system_a_saem_sources.sha256  pinned external SAEM inputs
slurm/system_a_saem_anchor.sbatch    all-115 anchor-only SAEM fit
slurm/system_a_anchor_banks.sbatch   12-patient bank array
slurm/system_a_candidate_banks.sbatch component-target bank array
slurm/system_a_message_validation.sbatch ODE-free plan/assessment
slurm/system_a_pure_anchor_banks.sbatch final pure-SAEM reference array
slurm/system_a_pure_candidate_banks.sbatch final validation-target array
slurm/system_a_pure_message_validation.sbatch final ODE-free plan/assessment
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

Generated logs and large RDS banks are excluded from Git. Compact result
summaries and review evidence may be committed when needed for audit.

## System A prerequisites

The System A pilot is not standalone. It requires the local certified
`modular_bayes` target mirror, numerical-audit certificate, and clean pinned
worktree at commit `a3e0367b06c82ad4d07280c99f10d4f9bac69978`. It does not
read comparator posterior draws. The adapter checks commit identity, source
hashes, target fingerprint, and the numerical-audit certificate before loading
the sealed VODE-BDF callbacks. The Slurm launchers contain Oxford/user-specific
absolute paths and must be adapted elsewhere.
