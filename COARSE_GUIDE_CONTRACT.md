# Coarse forward guide: a gated route towards reusable patient modules

The goal is modular Bayesian inference across population models. This first
test asks only whether a cheap forward calculation could accelerate exact
patient refreshes at the CURRENT shared parameters. It cannot establish
marginal-message accuracy, evidence, full shared-parameter uncertainty or useful
cross-model reuse. No main-project or System B source is changed.

## Exact correction

At fixed globals c, write the exact patient factor F(x;c)=L(x;psi)g(x;eta).
The cheap target uses the SAME population density and a strictly positive guide
likelihood: F_tilde(x;c)=L_tilde(x;psi)g(x;eta). Repeating one reversible cheap MH
kernel a fixed k times gives a reversible macro proposal. Its endpoint correction
is exp(min(0, [log F-log F_tilde](x')-[log F-log F_tilde](x))). Rejected proposals
leave both exact and approximate active caches unchanged. No stationarity of
the finite inner subchain is required. No Gaussian population assumption is
required; the exact g is evaluated explicitly. Changing globals requires caches
to be refreshed consistently. The guide is fixed during a retained epoch.

This is the established surrogate-transition construction, not a new identity:
Lykkegaard et al., section 2.2.1,
https://epubs.siam.org/doi/10.1137/22M1476770 . Its accuracy and mixing depend on
the guide. It does not make shared-dynamic updates or all later models ODE-free.

## Stage 1: measure cheapness before committing to chains

- Three already-audited ODE patients: 3, 74 and 122, all eight local coordinates.
- Ninety-six predeclared design states per patient: first 96 independent
  validation states from `outputs/residual_ceiling_122250/design_and_states.rds`,
  SHA256 `f350728c91de7b568ff50e44ece79408d9c7682c534bde473b5716d8c1125c4c`.
  These are Gaussian proposal-design states, NOT patient conditional draws.
- Three shared contexts: the saved oracle centre, +one whitened unit in dynamic
  direction 1, and +one whitened unit in direction 2. Whitened directions come
  from `outputs/system_a_12_oracle_rescue_v2/frozen_map.rds`, hash
  `ec4b29a4efc7162da43f2295c0530d1f7a02e72e1c785df1cd0a7e3b6b590a78`.
  Population densities are irrelevant to this pointwise numerical check.
- Three fixed VODE-BDF guide profiles: (rtol,atol)=(1e-4,1e-6), (1e-2,1e-4),
  and (1e-1,1e-3). Only tolerances change in a private patient-data copy.
  Equilibria, observation times, censoring and heteroscedastic observation
  likelihood are supplied by the existing pinned implementation.
- The original adapter remains the sealed target, with rtol=1e-6 and atol=1e-8.
  Guide likelihood failures are floored at log likelihood -1e6 for proposal
  support only. Approximate values never enter an exact cache.
- Two fresh evaluations per state and profile, including the exact profile.
  Rotate execution order to reduce systematic warm-cache timing bias. Check
  deterministic repeat equality. Count CPU time, elapsed time, forward attempts,
  actual ODE integrations and failures, including non-finite exact states.
- At most 1,728 exact and 5,184 approximate prediction calls. No new bank,
  gradient training, fitted density or neural network.

The preflight asks for a single profile with median CPU speedup >=2 across the
nine patient/context cases, speedup >=1.5 in every case, >=90% finite exact
design states, no approximate failure on an exact-finite state, and 95th
percentile absolute DIFFERENCE in log-likelihood residual <=1 in every case.
Residual differences use consecutive disjoint design pairs and are reported
as numerical-guide diagnostics, not stationary sampler acceptance. These are
predeclared exploratory go/no-go gates, not universal efficiency theorems.
If several profiles qualify, take the fastest median profile and freeze it
before fresh sampling. If none qualifies, stop this coarse-tolerance route
before constructing a larger sampling experiment.

## Stage 2, only after a qualifying guide exists

Compare a full-dimensional exact local MH baseline with the same proposal run
in cheap reversible subchains plus exact endpoint correction. Freeze separate
seeds, warmup, proposal tuning and budgets before running. Use multiple chains
and patient population-score summaries, not acceptance alone. Charge all
guide calls, failed proposals, setup and elapsed time. A positive result needs
conditional posterior agreement, convergence and a material ESS/CPU gain.

Even that would only establish an accelerator. Progress towards MODULAR Bayes
then requires reuse across meaningfully different population structures, with
scientific densities and covariate effects explicit, full joint uncertainty,
and cumulative new ODE cost versus separately fitting the structures. No
promotion or 115-patient run follows automatically from either pilot stage.
