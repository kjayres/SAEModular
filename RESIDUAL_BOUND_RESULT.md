# Residual-collapse bound feasibility: no-go for the supplied construction

Job `122247` completed successfully on 2026-09-05. This is a structural audit
of one proposed construction, not a failed run of a complete residual sampler.
The [frozen specification](RESIDUAL_BOUND_FEASIBILITY.md) states the scope and
analytic arguments. [Compact evidence](results/residual_bound_probe) is committed.

## What the check established

- **An everywhere-positive Gaussian lower component is incompatible with the
  sealed support.** Four finite-coordinate witnesses (two at each frozen
  context) had finite population log density but exact log likelihood `-Inf`
  because their equilibrium flux was negative. A strictly positive B cannot
  satisfy B <= L there.
- **The proposed finite global affine-forward error assumption fails.** Baseline
  CD4 grows exponentially in log lambda along a valid parameter ray. Its error
  against the anchor tangent has the analytic form
  `C*lambda_A*(exp(d)-1-d)`, and is unbounded. The sealed baseline predictions
  agreed with this expression at all twelve tested positive displacements.
  At d=1 the error was about 2.28 and 2.50 fixed-anchor CD4 standard deviations;
  at d=4 it was about 157 and 173. These noise units are diagnostic only: actual
  System A CD4 error SD changes with its predicted mean.
- **The pasted observation-model assumptions also do not match System A.**
  Patient 3 includes a censored viral observation and prediction-dependent CD4
  error SD, rather than an uncensored fixed-covariance Gaussian likelihood.

These analytic obstructions concern the global recipe. They do not show that a
locally useful approximation, a different log-likelihood bound, or a richer
support-aware residual construction is impossible.
The unbounded-error argument concerns the mathematical forward model in the
current transformed coordinates. Restricting to a finite region or to numerical
solver-success states is a different premise requiring its own support and
error control; the finite numerical checks are not themselves a global proof.

## A limited constructive finding

Simple intersections of two linear half-spaces guarantee a valid mathematical
equilibrium. Their probability under the current diagonal Gaussian population
is cheap and analytic. For the predeclared alpha-based family and grid:

| Frozen context | Best k | Region population mass | Estimated full valid-equilibrium population mass |
|---|---:|---:|---:|
| Comparator chain 1 | 1/128 | 0.614% | 98.682% |
| Comparator chain 3 | 1/256 | 0.812% | 99.487% |

The last column uses 4096 independent population draws per context, not patient
posterior draws. No sampled point inside these regions violated the equilibrium
condition. Four points from each selected region also passed the sealed solver.

The small region masses indicate conservatism of this particular family, not
that most population states are invalid. Other regions could be better. Neither
population coverage nor sampled solver success certifies a useful likelihood
lower bound, a collapsed posterior fraction b/M, or numerical solver success
throughout a region. No such claim is made.

## Cost and decision

The audit used **26 sealed prediction calls: 22 actual ODE integrations and four
pre-integration equilibrium rejections**. There were no integration failures.
The measured probe phase took about 4.25 seconds; the complete one-core Slurm job,
including startup/loading, took one minute. No SAEM, patient MCMC bank, fitted
transport, SMC, neural network, or retained residual chain was built.

**Decision:** do not implement the full residual-collapse sampler using the
supplied Gaussian/affine construction. A useful lower likelihood component,
including support/tail handling and an affordable integration/sampling interface,
is still missing. The support repair alone is insufficient. Park sampler
implementation; any continuation should first solve that lower-bound problem.

This does not falsify the residual-augmentation identity or all possible bound
families. The separate exact-recycling comparison remains a different experiment.
