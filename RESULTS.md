# Experiment status

All ODE experiments below ran as Oxford Statistics Slurm jobs. Generated
patient-level artifacts remain outside Git.

## Nonlinear toy

Final validation job: `121937` (7 minutes 7 seconds, one CPU).

The exact affine chain agreed with the high-order quadrature posterior and
passed the frozen rerun gate.

| Diagnostic | Fixed states | Affine transport |
|---|---:|---:|
| Global acceptance | 0.441 | 0.383 |
| Total initial-positive-sequence ESS | 1374.5 | 3678.3 |
| ESS per million total patient-likelihood calls | 596.6 | 1596.4 |
| Posterior-scaled ESJD | 0.293 | 0.521 |
| Split R-hat | 1.0044 | 1.0004 |

Observed affine/fixed ratios were 2.676 for ESS per call and 1.780 for ESJD.
The affine chain's posterior mean was 0.0084 reference standard deviations from
the quadrature mean. The map-only cost amortized after approximately 640
effective shared-parameter samples when a SAEM fit was assumed already
available; counting SAEM and map construction from scratch gave approximately
945 samples. The selected affine proposal was at the largest tested scale, so
these are validated observed results, not a claim of fully optimized tuning.

## System A: initial 12-patient oracle pilot

Initial job: `121936` (4 minutes 50 seconds, 12 CPUs). Source snapshot:
`761f2ac`. Certified target fingerprint:
`123fc00a2c0151215f9c550613819d4563a02c5bbe789eee5f2aee2c20844925`.

This first pilot failed both its map-quality and mechanism gates.

| Proposal radius | Fixed acceptance | Affine acceptance | Affine/fixed ESJD per ODE |
|---:|---:|---:|---:|
| 0.5 | 0.651 | 0.545 | 0.837 |
| 1.0 | 0.410 | 0.221 | 0.538 |

At radius 1, affine total work variance ranged from 4.264 to 6.020 and forward
D2 from 3.349 to 6.287 across the four signed axial moves. There were no replay
solver failures. However, minimum conditional-bank coordinate ESS was only 4.4,
so the 36 Cholesky quantities and their global sensitivities were not estimated
credibly. The original central-difference map also failed to reproduce the
sampled endpoint moments under curvature or Monte Carlo error.

Decision: do not run 115 patients or retained System A chains. Run one bounded
12-patient rescue with longer independent fit/held-out banks, interpolation that
passes exactly through the axial training moments, and a conditional-mean
translation comparator. If adequately mixed banks still fail to beat fixed
states, stop this dynamic-shared-parameter transport route for System A.

This System A result is not evidence yet for a SAEM handoff, population-level
posterior uncertainty, bridge sampling, retained-chain correctness, or
patient-count scaling.

## System A: bounded 12-patient rescue

Rescue job: `121938` (27 minutes 10 seconds, 12 CPUs). Source snapshot:
`cb34c71`. The job completed with exit code zero and all target and exact-kernel
tests passing.

The predeclared mechanism gate did **not** pass. The conditional banks did pass
their quality gate, so this is more informative than the initial pilot.

| Proposal radius | Fixed acceptance | Mean-only acceptance | Affine acceptance | Mean/fixed ESJD per ODE | Affine/fixed ESJD per ODE |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.648 | 0.756 | 0.697 | 1.165 | 1.075 |
| 1.0 | 0.408 | 0.557 | 0.452 | 1.366 | 1.107 |

The continuation threshold was fixed at 1.5 times the fixed-state ESJD per ODE.
Mean-only transport was the best method but reached only 1.366; full affine
transport reached 1.107. The quality and numerical diagnostics were otherwise
sound:

- minimum coordinate ESS across map and held-out banks: 56.4;
- maximum independent-centre whitened mean discrepancy: 0.431;
- maximum held-out covariance discrepancy in whitened Frobenius norm: 0.828;
- maximum axial training-moment reconstruction error: `3.33e-16`;
- maximum cross-transport log-normalizer disagreement at radius 1: 0.117;
- no mean-only or affine endpoint solver failures.

The gain was directional. At radius 1, mean-only transport was essentially
neutral for the first whitened shared-parameter axis (0.990 times fixed) and
2.179 times fixed for the second axis. Mean transport reduced patient work
variance for 11 of 12 patients, so that second-axis signal was not created by a
single patient. Full covariance transport was consistently worse than moving
the conditional mean alone. The useful second-axis diagnostics were stable
across half-banks; a reverse-D2 half-bank discrepancy of 0.802 in the
first-axis-minus move is an additional reason not to over-interpret the full
average.

Map and replay-bank construction cost 936,072 exact patient ODE calls, while
the paired endpoint replay used 864,000 proposed ODE calls. This startup cost is
far too large to justify a 1.366-fold endpoint gain unless the mean map can be
obtained nearly for free from a preceding SAEM computation. Even excluding the
held-out validation bank and anchor cost, the observed gain would require about
190,000 one-standard-deviation mean co-move attempts to amortise the present
map-fitting calls. A purpose-built mean-only fit could be cheaper, but was not
tested here.

Decision: do not run the 115-patient affine experiment or a retained System A
chain. Freeze the full-affine branch. The next small System A test should target
the hierarchy-conditioned local patient proposal at the current shared state,
which directly addresses the audited 3.62% local-update acceptance. Mean-only
global transport remains a possible later component for the second shared
direction if a genuine SAEM handoff supplies its map cheaply.

The rescue is still an oracle-anchor endpoint replay, not a SAEM implementation
or a posterior sampler. It provides evidence for a limited conditional-mean
effect, but not for fast population-level uncertainty or patient-count-invariant
efficiency.

## System A: corrected all-patient SAEM anchors

Corrected jobs `122128` (seed 8) and `122129` (seed 29) both used all 115
patients, the 300 + 100 iteration schedule, and one SAEM chain. They completed
in 12 hours 12 minutes and 12 hours 1 minute, respectively. The initialization
overlay worked as intended: the eight active random-effect variances in the
saved iteration-zero parameter vector exactly equal the published starts, the
four initialization-only dummy variances are `1e-6`, and the corresponding
fitted structural variances remain zero. Both jobs captured zero R warnings.

The fitted endpoints are nevertheless incompatible.

| Quantity | Seed 8 | Seed 29 |
|---|---:|---:|
| `u_eta_rti` | 31.259 | 14.087 |
| `u_eta_pi` | 7.487 | 39.835 |
| `beta_nelf` | -11.421 | -28.064 |
| `u_log_mu_t` | -5.845 | -7.884 |
| `omega_mu_t` variance | 2.643 | 4.483 |
| Invalid-equilibrium proposals | 18,782 | 20,061 |

The implied untreated/treated protease-inhibitor efficacies are approximately
`0.99944 / 0.01919` for seed 8 and `1 / 0.999992` for seed 29. RTI efficacy is
also effectively one in both fits. Several random-effect variances either
collapse or differ materially between seeds. All saved final patient states
are finite and have analytically valid untreated equilibria, but the corrected
fits have not yet been converted through the sealed VODE-BDF handoff ledger.
No likelihood objective was requested, so there is no defensible basis for
selecting one endpoint after inspection.

Decision: the software initialization defect is fixed, but neither result is a
credible sole SAEM anchor. Do not average them or choose the more convenient
one. They may be retained as predeclared stress branches in a bounded
multi-anchor experiment; otherwise the model/fit requires stronger
identification before a sole-anchor claim.

## System A: defensive-reference message calibration

Reference-bank array `122130`, ODE-free plan `122142`, component-bank array
`122144`, and recorded assessment `122170` used the old production SAEM anchor
on the 12 audited patients. The reference was

```text
h_i = 0.5 g_SAEM,i + 0.5 g_priorcentral,i.
```

The assessment status is `inconclusive_calibration_or_mixing_failure`. The
single defensive-reference chain is not usable here:

- mean reference-bank sampling acceptance was 9.39%, with 24 of 48 chains
  below 2%;
- the worst coordinate ESS was 2.80, and no chain met the predeclared ESS 400
  requirement;
- the prior-component raw message had minimum relative weight ESS 0.001 and a
  maximum normalized weight of 1;
- the prior-component pooled bridge failed for patients 3, 6, 55, and 74;
- its cohort forward chain range was 1774.47;
- the SAEM-component cohort forward estimate was 1.757 versus reverse 8.318,
  also failing overlap, replication, and agreement gates.

The separated component-target pCN chains were better computationally (about
42% mean acceptance for the SAEM component), but their minimum per-chain
coordinate ESS was still only 6.44 from 1,000 retained draws. The component
closure identity is algebraic and does not independently validate sampling.

Decision: stop the SAEM-plus-prior mixture as a single-bank strategy and do not
scale it to 115 patients. This experiment did not evaluate genuine nearby
population endpoints from a pure SAEM bank, so it does not by itself reject
that narrower idea. The only justified continuation is one bounded
12-patient, pure-SAEM, fixed-shared-parameter test at the corrected SAEM
branches, with a strict post-SAEM ODE budget and independent candidate-bank
validation. Failure with adequately mixed banks stops the single-anchor
proposal; success would justify, but not itself complete, a 115-patient test.
