# SAEModular

This research prototype tests one idea: use an all-115-patient SAEM fit as a
hierarchy-aware anchor, then turn independent exact-target patient-conditional
MCMC banks into reusable messages for fast modular population inference. The
first reference broadens the fitted region toward one explicitly declared
prior-centred region; it does not claim to cover the prior's full support:

```text
h_i(x) = 0.5 g_SAEM,i(x) + 0.5 g_priorcentral,i(x).
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

The first System A implementation is narrower than this generic identity. Its
two reference components are eight-dimensional diagonal Gaussians whose means
include the observed Nelfinavir indicator, but `h_i` itself is a non-Gaussian
mixture. An independence proposal from `h_i` can cross between its separated
components; because the proposal and reference cancel, its exact MH correction
is only the sealed likelihood ratio. It does **not** require
candidate message models `g_(i,M)` to be diagonal Gaussian: those models enter
later only through their normalized log-density evaluations. An arbitrary
non-Gaussian anchor would require a different exact-target conditional MCMC
kernel or direct reference sampler, but not a different message estimator.
Likewise, bidirectional bridge
validation for an arbitrary candidate requires a valid candidate-target bank
sampler; the current bank builder covers only System A's diagonal-Gaussian
population family.

The repository now exposes the distribution-free log-density-ratio estimator,
but it does not yet contain a catalogue of candidate `g_(i,M)` models or the
outer sampler over `(M, eta)`. The first message experiment must supply those
normalized candidate log densities explicitly and keep their covariate data
and model priors under provenance control.

The first falsification sequence is:

1. fit the genuine all-patient SAEM anchor;
2. freeze all four shared parameters at that anchor;
3. build independent exact-target VODE-BDF conditional MCMC banks for the 12
   audited patients under the equal defensive reference;
4. calculate raw `g_SAEM / h` and `g_priorcentral / h` messages without ODE
   calls, and audit `0.5*g_SAEM/h + 0.5*g_priorcentral/h = 1` pointwise;
5. validate both separated components using separately simulated candidate
   banks and candidate-augmented bidirectional bridge estimates (the bridge
   and raw estimates share the held-out reference draws, so their agreement is
   not an independent replicate check);
6. only after that calibration passes, test treatment-effect and non-boundary
   scale directions;
7. proceed to all 115 patients only if overlap, replication, and accuracy pass.

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
experiments/toy_decay.R        nonlinear toy validation
experiments/system_a_12.R      12-patient oracle endpoint pilot
experiments/system_a_anchor_bank_patient.R  one-patient exact-target MCMC worker
config/system_a_saem_sources.sha256  pinned external SAEM inputs
slurm/system_a_saem_anchor.sbatch    all-115 anchor-only SAEM fit
slurm/system_a_anchor_banks.sbatch   12-patient bank array
slurm/system_a_candidate_banks.sbatch component-target bank array
slurm/system_a_message_validation.sbatch ODE-free plan/assessment
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

The System A pilot is not standalone. It requires the local certified
`modular_bayes` target mirror, numerical-audit certificate, and clean pinned
worktree at commit `a3e0367b06c82ad4d07280c99f10d4f9bac69978`. It does not
read comparator posterior draws. The adapter checks commit identity, source
hashes, target fingerprint, and the numerical-audit certificate before loading
the sealed VODE-BDF callbacks. The Slurm launchers contain Oxford/user-specific
absolute paths and must be adapted elsewhere.
