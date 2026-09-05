# Conditional sampling comparison, frozen after numerical preflight

Preflight 122252 qualifies profiles rtol=0.01 and 0.1. Following the original
selection rule, freeze the fastest median profile: VODE-BDF rtol=0.1,
atol=0.001, log guide floor=-1e6. Its median measured speedup was 3.099, minimum
2.691 and worst design-pair 95% absolute residual difference 0.935. Those are
not retained-chain results. The exact solver and target remain unchanged.

## Scope and target

Three patients (3,74,122), three fixed shared contexts (oracle centre and the
two +one-whitened-direction contexts from preflight), four independent chain
seeds per method/context/patient. These are 72 short conditional chains, NOT
72 patients or a joint population posterior. Population parameters are fixed
at the saved oracle anchor. All eight patient coordinates move. Both targets
are the original likelihood times the original population density, including
the treatment/covariate mean effect. System B is untouched.

The methods are a locally preconditioned random-walk MH baseline and the SAME
proposal run as eight cheap reversible MH steps, followed by one exact endpoint
correction. This isolates the benefit of cheaper inner exploration. It is not
a comparison against a deliberately weak untuned pCN baseline. The Gaussian
proposal is a computational choice, not an assumption that every future
population distribution is Gaussian; g is evaluated explicitly in both targets.

## Fixed run sizes and starts

- Baseline: 6,000 transitions, discard first 500.
- Guide macro: 1,500 transitions, discard first 125 (1,000 cheap inner steps).
- Four independent seeds per arm; matched initial states across methods only.
- Proposal covariance: saved centre covariance from the old conditional-map
  artifact, with initial multiplier 2.38/sqrt(8). No transport or covariance
  derivative is used. Scalar step size adapts ONLY during each warmup, targeting
  inner acceptance 0.234, and is frozen for retained draws. This is an oracle
  local-efficiency test; historical covariance estimation is excluded from both
  costs and cannot be claimed as a free deployable preparation step.
- Starts use unused source validation rows 193,257,321,385 per patient, with
  up to 16 sequential fixed candidates for each chain if an exact state is
  invalid at that shared context. All attempted initial evaluations are charged.
- At most 271,728 exact calls and 432,036 cheap calls across the actual chains,
  excluding a separately labelled implementation smoke. Actual calls can be
  lower if the inner endpoint is unchanged. All failures count. Fixed transition
  counts avoid stopping a chain based on its observed convergence or CPU cost.
- Up to four Slurm workers. Save every retained original patient vector, exact
  log likelihood, tuning and per-category counters/times. Checkpoint every 500
  outer transitions. No new density fitting or surface construction.

These counts produce comparable, not identical, anticipated CPU budgets based
on the preflight speed ratio. The PRIMARY comparison uses actual total CPU time,
not raw iteration counts or exact-call counts alone. Include warmup, initial
evaluations, bookkeeping and setup. Charge all 40.182 preflight CPU seconds
(Slurm accounting for job 122252, including tested but unselected profiles) to
the guide arm. Split common preparation across arms. Retain full cost ledgers.

## Diagnostics and decision

For each patient/context, compare four-chain distributions for the eight local
coordinates, eight normal-population variance scores z_j^2-1, and exact log L.
Coordinate autocorrelation is identical to the corresponding scaled population
mean-score autocorrelation at fixed globals; the treatment coefficient's score
is proportional to coordinate eight for treated patients and zero otherwise.
No acceptance-only or raw-coordinate-only mixing claim.

Require R-hat<=1.01, bulk ESS>=400 and tail ESS>=200 for every diagnostic.
Compare means, SDs and 5%,50%,95% quantiles with combined MCSE and a Bonferroni
95% family tolerance across all contexts/quantities/statistics. Report score ESS
per CPU, wall time and exact call, plus IACT and lag-one persistence.
Performance go: median score bulk ESS/CPU ratio >=1.5 overall, median ratio >=1
in every patient/context, and no score ratio <0.8. A go also requires all
convergence and agreement checks. Nonconvergence makes the posterior comparison
inconclusive, not evidence of target bias. Do not automatically lengthen failures.

If this passes, the NEXT scientific test is cumulative cost and posterior
accuracy across different population structures, with shared uncertainty restored.
This experiment does not establish modular messages, Bayes factors, 115-patient
performance or full dynamic-shared uncertainty. A fast patient kernel is only
an enabling component for that overarching modular-Bayes goal.
