# Residual-collapse prerequisite: System A

Frozen before this probe, 2026-09-05. This is a structural feasibility check,
not another sampler, fitted surrogate, SAEM run, or performance demonstration.
The residual-augmentation identity is valid conditional on its assumptions.
We test the proposed concrete Gaussian/affine construction, not every possible
lower bound. The upper envelope is optional when using the exact-MCMC fallback;
a useful lower component and its integration/sampling interface are not optional.

## Predeclared expectation and scope

We expect the everywhere-positive Gaussian lower component and the finite-global-
error affine-forward recipe to be inapplicable. The reasons below are analytic;
the numerical probe checks them against the actual adapter. A support-aware
repair is investigated separately, without calling it a likelihood certificate.

Use ODE patient 3 at both previously frozen comparator contexts. Read the existing
probe design with SHA256 51fc4d28614e42563f9b0f2f460a77f9b630e3a96f5920bc17933d4189ce3fd8.
No new fit or conditional bank. At most 40 sealed prediction calls; 4096 cheap
population/equilibrium checks per context. All numerical work runs on Slurm.
Report all attempted calls and failures. Do not modify the running recycling
experiment, System B, or the sealed numerical target.

## Two obstructions to the pasted construction

Write x for the eight canonical transformed patient coordinates. With all other
coordinates fixed, the equilibrium has T0 independent of lambda and flux
lambda - mu_t*T0. Values lambda below mu_t*T0 form an open invalid region with
positive Gaussian population probability. The sealed likelihood there is zero.
Any nonzero, everywhere-positive exponential quadratic B violates B <= L there.
Re-evaluate lambda at 0.5 and 0.9 times the threshold to record concrete witnesses.

Patient 3 has baseline CD4 observations. Put

    C = (1-pi)/(alpha_l+mu_l)
        + (alpha_l+pi*mu_l)/(mu_a*(alpha_l+mu_l)) > 0.
    CD4(0) = T0 + C*(exp(x_lambda) - mu_t*T0).

Consequently any affine forward map in x has unbounded absolute error along
x_lambda -> infinity, even within the valid equilibrium region. With fixed
positive-definite observation precision, no finite global error constant eta
of the proposed affine-forward kind exists. For the anchor tangent, the exact
baseline error at displacement d is C*lambda_A*(exp(d)-1-d). Check d in
{0.1,0.25,0.5,1,2,4} against the sealed baseline predictions. ODE failures remain
recorded, not treated as successful likelihood evaluations.

The observation model also includes censored viral measurements and CD4 standard
deviation proportional to its predicted mean. It is not the fixed-covariance,
uncensored Gaussian model assumed by the proposed closed-form bound derivation.
The fixed anchor noise scale in the tail CSV is ONLY an illustrative scaling,
never a replacement observation model or a certified error constant.

## A constructive support repair to test

The support obstruction alone need not destroy cheap Gaussian integration.
Define

    a = x_lambda + x_p - x_mu_t - x_mu_a + log_gamma - log(mu_v).
    D_k = { x_alpha >= log_mu_l + log(k), a >= log(1+1/k) }, k > 0.

For alpha_l >= k*mu_l,

    (alpha_l+mu_l)/(alpha_l+pi*mu_l) <= 1+mu_l/alpha_l <= 1+1/k.

Therefore D_k guarantees a valid mathematical equilibrium (nonnegative flux).
It is the intersection of two linear half-spaces in x, not a nonlinear
truncation. Under the current diagonal Gaussian population, a and x_alpha are
independent normals, so P_g(D_k) is a product of two normal tail probabilities.
Under an arbitrary nonsingular Gaussian density the corresponding probability
is a bivariate Gaussian probability. Thus multiplying a Gaussian-shaped
component by this support indicator need not require an eight-dimensional
normalization. This statement does NOT extend automatically to arbitrary
non-Gaussian population distributions.

Predeclare k=2^j for j=-8,...,8. Report analytic population mass and Monte Carlo
membership checks at each context; no likelihood calls are used to select k.
Check four population draws in the largest-mass region with the sealed solver.
The equilibrium proof is analytic; agreement at sampled points is not a global
ODE-solver-success or numerical-likelihood certificate. Overflow and solver
failure regions remain an additional issue for the sealed numerical target.

Most importantly, D_k does not provide the amplitude/shape of B <= L *inside*
the region. Its population mass is NOT a collapsed posterior fraction, b/M, or
a speedup estimate. No integrated residual-gap kappa is claimed without such B.

## Decision

An invalid-state witness rules out the unmodified positive Gaussian lower bound.
The analytic exponential baseline tail rules out the proposed global affine-
forward error assumption in these coordinates. Neither rules out all residual
collapse methods. Cheap support regions, if confirmed, retain a possible avenue
but do not authorize a full sampler. A further step would first need an explicit
useful lower likelihood component, correct tail/numerical support, and affordable
integrals and draws under the population query families of interest.

Keep this probe and its compact results; do not launch a hierarchy or add a
certification framework on the basis of this audit.
