# Local residual component: bounded rejection screen

This is a diagnostic, not a residual sampler or a likelihood certificate.
The global Gaussian construction remains retired. The new question is whether
a specified compact component could collapse a useful fraction of one patient's
conditional target. No new SAEM fit or patient MCMC bank is required.

## A bank-free upper confidence bound

At a fixed global context let F be the complete unnormalised patient factor,
including the original population density and any coordinate Jacobian, and
M = integral F. For a fixed normalised candidate q, a valid lower component
c q <= F must have c <= c_probe = min_j F(x_j)/q(x_j), at interior probe points.
An exact zero at such a point rules out a positive pointwise lower multiplier.

Independently of fitting and probing, draw N independent states from any fixed,
normalised density h. For each predeclared cap C > 1, form

    Z = min{ F(X)/(h(X) c_probe), C },  X ~ h.

Then E[Z] <= M/c_probe. Hoeffding's inequality gives the lower confidence bound

    L = mean(Z) - C sqrt{ log(1/alpha)/(2N) }.

When L > 0, every valid component therefore satisfies, with probability at least
1-alpha, c/M <= min(1, 1/L). Otherwise the upper bound is the vacuous value one.
The validation density h need not equal q: one shared validation set can assess
several regions. Clipping makes the concentration argument finite without an
upper likelihood bound or a finite variance assumption on F/h.

This is a conservative Monte Carlo diagnostic with an explicit error probability,
not a deterministic certificate. It is not an importance-sampling replacement
for the posterior. In particular a high bound proves neither coverage nor a
useful lower multiplier, and a low bound rejects only the specified q. Floating-
point numerical evaluation is treated consistently with the sealed target, not
as interval-certified arithmetic or exact continuous-ODE integration.

## Frozen first screen

- Patients 3, 74 and 122: three already-audited ODE patients, including both
  treatment groups. They are examples, not a representative cohort sample.
- One fixed oracle context: the saved successful conditional-map fitting run
  `outputs/system_a_12_oracle_rescue_v2/frozen_map.rds`, SHA256
  `ec4b29a4efc7162da43f2295c0530d1f7a02e72e1c785df1cd0a7e3b6b590a78`.
  Only its numerical centre/covariance coefficients are reused. No transport
  kernel or archived implementation is restored; no old draws are presumed
  stationary for a different conditioning context.
- Keep all eight local coordinates. Replace log(lambda) by
  u = log(lambda/lambda_crit - 1), with the other seven unchanged. Preserve the
  original population density and multiply by d log(lambda)/du = logistic(u).
  No Gaussian population distribution is assigned to u and no population
  density is renormalised over feasible equilibria.
- Transform the saved conditional mean and covariance locally using the exact
  coordinate Jacobian. These moments define a candidate, not a certified fit.
  An infeasible saved mean is reported rather than silently repaired.
- q is this Gaussian truncated to its own Mahalanobis ball, with Gaussian
  probabilities 0.5, 0.8 and 0.95. Radii use chi-squared probabilities in eight
  dimensions, not one-dimensional standard-deviation coverage.
- h is the untruncated same Gaussian in u coordinates. Per patient, 512 fresh
  independent h draws serve all three regions. No MCMC convergence premise.
- Each q uses 64 independent q draws, the centre, and 16 axial probes at
  0.999 times the region radius. Probe and validation streams are independent.
- Caps 2, 10 and 50. Allocate alpha = 0.05/(3 patients * 3 regions * 3 caps),
  giving simultaneous 95% coverage by a union bound; dependence across the
  different reported bounds does not invalidate that argument.
- At most 2,500 sealed prediction calls including all probes and validation
  (planned 2,265). Count actual ODE integrations, pre-solver rejections, solver
  failures and time separately. Save numerical states and likelihood outputs
  needed to audit the calculation, without generating new MCMC traces.

## Decisions

The practical screen is whether the upper confidence bound rules out 80% collapse.
That is a predeclared aspiration for substantial work reduction, not a universal
efficiency theorem: active-state mixing and certification costs remain untested.
If all three regions fail for a patient, retire this fitted local Gaussian family
for that patient. Do not infer failure of all compact components or all patients.
Any survivor is merely unresolved and does not authorise a residual sampler or
a 115-patient expansion. A later continuation would require one affordable,
valid lower certificate and a materially different justification, not automatic
bank extension or tuning to this validation set.

The independent shared-recycling retained-chain experiment is reported separately.
