# Local compact residual components: negative bounded screen

Job `122250` completed successfully on 2026-09-05. The experiment used all eight
local coordinates for three audited patients at one fixed oracle context. It
tested the new equilibrium-excess coordinate with the original population
density and the correct Jacobian, not a Gaussian replacement population model.
The [contract and derivation](RESIDUAL_CEILING_CONTRACT.md) were fixed before
submission. Thirty implementation assertions passed before the experiment.

## Result

The numbers below are **simultaneous 95% upper confidence bounds on the largest
possible collapsed posterior fraction** for each specified component. They are
not attained fractions, posterior coverage estimates or lower-bound certificates.
The row label describes the mass of the *fitted Gaussian* inside its truncation
ball, not the mass of the true patient posterior.

| Gaussian mass in component region | Patient 3 | Patient 74 | Patient 122 |
|---|---:|---:|---:|
| 50% | 22.82% | 29.32% | 14.82% |
| 80% | 17.93% | 38.08% | 12.86% |
| 95% | 3.63% | 17.79% | 5.52% |

All nine components fail the predeclared aspiration of 80% collapse. Even taking
the most favourable of the three radii for each patient leaves upper ceilings of
22.82%, 38.08% and 14.82%. Thus these components cannot marginalise most patient
updates, even before charging the cost of obtaining a valid lower certificate.

The negative outcome is not an MCMC mixing failure or a consequence of having
the wrong conditional bank: no patient chains were used for validation. Nor is
it the earlier invalid-equilibrium Gaussian obstruction: the new chart was used,
and all 2,265 endpoint integrations succeeded. The numerical geometry of these
particular components still permits only small pointwise lower multipliers.

## What was measured, and why the conclusion is conservative

The saved conditional-map means and covariances supplied candidate geometry only.
They were not taken as an exact posterior density. Each region's multiplier was
bounded above by 81 independent/design probe locations. A separate 512 IID draws
per patient from a normalised, untruncated Gaussian h supplied a lower confidence
bound on the normalising integral, using bounded weights and Hoeffding's inequality.
The same validation states served all three regions; a union bound over all
patients, regions and caps handles their dependence and the minimum over caps.

This differs from the pasted posterior-bank ceiling: it is a one-sided, bank-free
rejection screen. It can be more conservative, but any low confidence ceiling
still excludes a large valid component with the stated error probability. No
claim of exact interval arithmetic or a certified continuous-ODE target is made.

These results concern a locally transformed Gaussian geometry estimated in the
earlier oracle experiment. They do not rule out a different geometry, a different
component family, other global contexts or residual-collapse identities generally.
The three patients are not a representative sample from which to extrapolate a
full-cohort collapse fraction. A moderate reduction in active mass might help some
other algorithm; this experiment does not measure cost per effective sample.

## Cost and decision

- 2,265 sealed prediction calls and 2,265 actual ODE integrations; no failures.
- About 18.8 CPU seconds and 39.1 elapsed seconds for the diagnostic phase.
- Complete one-core Slurm job: 1m 13s, including tests and setup.
- No new SAEM, patient MCMC bank, certification engine or residual sampler.
- Historical oracle map construction is not charged to this incremental screen;
  it remains part of any honest end-to-end deployment comparison.

**Park this local Gaussian residual-collapse route. Do not build its sampler,
scale to 115 patients or automatically refit a sequence of richer components.**
The useful deliverable is the small rejection diagnostic: it prevented an
expensive certification/sampling programme for unsuitable candidate components.
The broader mathematical construction is unresolved, not disproved.

[Compact evidence](results/residual_ceiling) contains the cap-wise bounds,
decisions, full call ledger and provenance. Original-coordinate and chart states
and frozen candidate geometry are retained locally in
`outputs/residual_ceiling_122250/design_and_states.rds` for audit.
