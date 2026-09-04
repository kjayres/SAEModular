# Final pure-SAEM System A falsification contract

Frozen: 2026-09-04, before any pure-SAEM patient bank was simulated.

## Question and scope

This is the last single-anchor System A test. It asks whether a short,
pure-SAEM-conditioned patient MCMC refresh supports accurate raw population
messages for meaningful nearby changes while all shared parameters are fixed.
It does not test a 115-patient posterior, dynamic-shared-parameter uncertainty,
relative mass of the two SAEM modes, new population structures, or exact
retained inference.

Both corrected SAEM endpoints are mandatory, independent stress branches:

| Branch | Corrected run | Canonical-anchor SHA-256 |
|---|---|---|
| seed 8 | `system_a_corrected_full_seed8_20260903a` | `ce601409d84bfe5492667e64e688b8a81b34947b8f5372be6767ff2d117130ec` |
| seed 29 | `system_a_corrected_full_seed29_20260903a` | `19630a7f11f7508a53939e94a89f12e76f91360bb5cac0e3291722a4cc3c559d` |

The same twelve previously audited patients are used. No branch may be chosen,
averaged, shrunk, or extended after inspecting its new banks.

## Stage 0: sealed handoff

Every one of the 115 saved SAEM patient states must pass one evaluation through
the sealed VODE-BDF System A target. Both branches passed this prerequisite in
jobs `122172` (seed 8) and `122171` (seed 29): zero callback errors and zero
ordinary exact-target rejections.

## Frozen patient-bank budget

Every patient-target bank uses pure pCN with respect to its diagonal System A
population distribution:

- four independent, dispersed chains;
- 250 warm-up proposals per chain;
- 500 retained proposals per chain;
- thinning one;
- initial pCN beta 0.4, adapted only during warm-up toward acceptance 0.3;
- at most 18 predeclared initialization candidates per chain;
- at most 3,072 exact prediction calls or ODE integrations per patient-target.

Reference base seeds are `840003` and `2940003`. Candidate jobs reuse the same
branch base, and the worker adds its frozen `500000` candidate namespace;
reference jobs use the separate `1000000` namespace. The deployable reference-bank cap
is 36,864 calls per 12-patient branch. Scaling the same budget to 115 patients
would cap reference construction at 353,280 calls per branch.

Candidate banks exist only for validation. If both Stage-1 branches pass, their
joint cap is 294,912 calls. Including reference banks, the complete two-branch
test is capped at 368,640 calls, plus the 230 Stage-0 target evaluations. No
chain extension or retuning is permitted after results are seen.

## Predeclared endpoints

Each branch freezes its own `psi` and every unlisted population coordinate.
Four endpoints are evaluated:

1. `beta_nelf - exp(log_omega_eta_pi)`;
2. `beta_nelf + exp(log_omega_eta_pi)`;
3. `log_omega_lambda - log(1.25)`;
4. `log_omega_lambda + log(1.25)`.

Thus the treatment mean moves by one fitted patient-population SD and the
lambda population SD changes by a factor of 0.8 or 1.25. The latter remains
inside the Gaussian raw-weight second-moment boundary because 1.25 < sqrt(2).

## Stage 1: ODE-free pilot gate

Only reference chains 01-02 decide whether candidate simulation is allowed;
chains 03-04 remain held out. Every endpoint at both branches must pass:

- relative raw-weight ESS at least 0.20;
- maximum normalized weight at most 0.05;
- within-chain half-message difference at most 0.25 log units;
- 12-patient chain-replicate range at most 0.50 log units;
- rank-normalized/folded split R-hat at most 1.01;
- bulk MCMC ESS at least 400 and tail ESS at least 200 for the raw-weight
  contribution and endpoint-relevant sufficient statistics;
- projected 115-patient log-message MCSE at most 0.50;
- no single pilot patient representing more than 50% of projected MC variance
  for a treatment endpoint, or 25% for a scale endpoint;
- every untreated patient's treatment-endpoint log weight, raw message, and
  bridge message equal to zero within `1e-8`, auditing the covariate identity;
- the exact frozen bank budget and artifact provenance checks.

For treatment endpoints, control-patient weights are exactly one, so the error
projection is `sqrt((35/4) * sum_treated(s_i^2))` over the four treated pilot
patients. For scale endpoints it preserves the full cohort's treatment mix:
`sqrt((80/8) * sum_control(s_i^2) + (35/4) * sum_treated(s_i^2))`.
Both projections assume the fixed panel is representative within its relevant
stratum; they are stress extrapolations, not empirical full-cohort standard
errors. A 25% concentration gate would be structurally attainable for four
treated patients only if their errors were exactly equal, so the treatment
gate is predeclared at 50% (at least two effective contributors) while the
12-patient scale gate remains 25%. Worst-patient sensitivity, the delta-method
log-bias proxy, and five-versus-ten-batch MCSE sensitivity are reported.

Failure of either branch stops the experiment before candidate-bank ODE work.

## Stage 2: independent validation gate

This stage runs only if both Stage-1 jobs succeed: both candidate arrays have a
joint `afterok` dependency on both plan jobs. It uses held-out reference chains
03-04 for forward raw importance, candidate chains 01-02 for reverse raw
importance, and reference chains 01-02 plus candidate chains 03-04 for bridge
estimates. Thus each side's bridge draws are disjoint from that side's raw
estimator, although the bridge-reference draws reuse the Stage-1 pilot bank.
Every bridge must
converge. It retains the Stage-1 0.20 relative-weight-ESS floor and the other
mixing, overlap, cost, and projected-error criteria, and requires patientwise
and cohort-level agreement among forward,
reverse, and bridge estimates. Cohort agreement has a minimum tolerance of
0.50 log units and otherwise uses twice the combined estimated MC uncertainty.

The experiment passes only if all four endpoints pass at both corrected SAEM
branches. A pass licenses only a 115-patient, fixed-psi validation of these four
local directions. Any failure closes the single-SAEM-anchor System A branch;
it must not trigger repairs with transport maps, defensive mixtures, flows, or
post-hoc smaller endpoints.
