# Final pure-SAEM System A result

Date: 2026-09-04

## Decision

The predeclared 12-patient falsification failed at Stage 1 under both corrected
SAEM anchors. Stop the single-anchor SAEM-message branch for System A. Do not
simulate candidate banks, run Stage-2 bridges, or scale this method to 115
patients.

This is a computational no-go under the frozen short post-SAEM budget. It does
not invalidate the importance identity. The observed raw-weight overlap was
adequate for all eight branch/endpoint combinations, but the pure-SAEM patient
chains did not produce reliable effective draws cheaply enough.

## Frozen run ledger

Source commit: `28920f234a90e66c8a396d2d6cbc102def096492`.

| Job | Role | Result |
|---|---|---|
| `122174` | complete source/synthetic test suite | passed, exit 0 |
| `122175` | seed-8 pure-SAEM reference array | all 12 tasks passed |
| `122176` | seed-29 pure-SAEM reference array | all 12 tasks passed |
| `122178` | seed-8 ODE-free Stage-1 gate | scientific fail, exit 3 |
| `122177` | seed-29 ODE-free Stage-1 gate | scientific fail, exit 3 |
| `122179`, `122180` | candidate arrays | dependency-blocked, then cancelled; no tasks ran |

Both canonical anchors had already passed all 115 sealed-target handoff checks
in jobs `122172` and `122171`. The exact design and stop rule are in
`FINAL_PURE_SAEM_CONTRACT.md`.

## Endpoint results

Every endpoint failed in both branches.

### Seed 8 branch

| Endpoint | min raw rESS | max split R-hat | min bulk ESS | min tail ESS | projected 115 MCSE | passed |
|---|---:|---:|---:|---:|---:|---|
| treatment minus | 0.350 | 1.065 | 23.30 | 22.38 | 0.896 | no |
| treatment plus | 0.271 | 1.093 | 23.30 | 22.38 | 0.938 | no |
| lambda scale ×0.8 | 0.812 | 1.175 | 7.29 | 15.33 | 0.340 | no |
| lambda scale ×1.25 | 0.809 | 1.175 | 7.29 | 15.33 | 0.282 | no |

### Seed 29 branch

| Endpoint | min raw rESS | max split R-hat | min bulk ESS | min tail ESS | projected 115 MCSE | passed |
|---|---:|---:|---:|---:|---:|---|
| treatment minus | 0.425 | 1.143 | 18.29 | 26.27 | 0.960 | no |
| treatment plus | 0.349 | 1.143 | 18.29 | 26.27 | 1.009 | no |
| lambda scale ×0.8 | 0.765 | 1.139 | 8.79 | 12.45 | 0.378 | no |
| lambda scale ×1.25 | 0.317 | 1.139 | 8.79 | 12.45 | 0.437 | no |

The gates were raw rESS at least 0.20, split R-hat at most 1.01, bulk ESS at
least 400, tail ESS at least 200, and projected full-cohort log-message MCSE at
most 0.50, together with split, replicate, concentration, cost, and identity
checks.

Raw overlap was not the immediate failure: all minimum rESS values exceeded
0.20 and all maximum normalized weights were below 0.05. The decisive common
failure was local MCMC persistence and between-chain disagreement. Relevant
bulk ESS was only 7.29--23.30 for seed 8 and 8.79--18.29 for seed 29 out of
1,000 pooled pilot draws. Maximum split R-hat was 1.175 and 1.143. Treatment
message uncertainty also projected to roughly 0.90--1.01 log units for 115
patients, about twice the 0.50 limit.

Other failures reinforce the stop. Seed 8 had maximum within-chain half-message
differences of 0.482, 0.921, 0.359, and 0.244 against a 0.25 gate; three cohort
chain ranges exceeded 0.50. Seed 29 had half-message differences of 0.593,
0.838, 0.285, and 0.522. Several projected-variance concentration gates also
failed. The untreated treatment-message identity and exact cost gates passed.

## Cost and interpretation

Each branch used 36,048 exact prediction calls for its 12 reference banks.
Successful ODE integrations were 26,796 for seed 8 and 26,897 for seed 29;
solver/nonfinite rejections were handled by the sealed target. Mean sampling
acceptance was respectable (0.419 and 0.427), showing again that acceptance did
not imply useful conditional mixing. Minimum single-chain coordinate or
squared-coordinate ESS was only 4.80 and 3.87.

Longer chains or a substantially stronger local sampler might eventually
estimate these messages, but that would be a new experiment and directly
weakens the proposed SAEM-cost advantage. Under the predeclared instruction to
stop when either branch cannot provide adequate mixing and cohort accuracy
within the short budget, no repair is licensed.

The conclusion is deliberately narrow: a short pure-SAEM pCN refresh did not
make the existing Stage-1/Stage-2 message identity into a fast System A modular
Bayesian method. Dynamic shared-parameter uncertainty, cross-mode posterior
mass, new structures, and a 115-patient posterior remain untested rather than
solved.
