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
