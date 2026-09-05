# Reduced exact hierarchical comparison

Frozen before retained-chain results, 2026-09-05. This is the separately
authorized follow-up to the shared-trajectory screen, whose numerical gate
failure remains recorded. Common-grid likelihoods now supply proposals only;
the retained target uses the original sealed System A numerical likelihood.

## Target and initialization

Use the twelve audited patients, in order:
3, 6, 16, 20, 55, 71, 74, 88, 105, 111, 117, 122.
The target is the original hyperprior, once, multiplied by these twelve patient
population densities and sealed likelihoods. All seventeen population and four
shared parameters remain uncertain. Do not temper the likelihood, rescale the
priors, or compare reduced results with the full-115 posterior as if it were the
same target.

Four independently seeded chains per method start from coherent comparator
draws: original chains 1--4, retained row 750. Each method receives the same four
starting states, restricted to the panel. These are valid initialization
candidates, not reduced-posterior draws. Freeze q to the equal mixture of the
treated/untreated population distributions at comparator snapshots 1 and 3.
The historical comparator computation is an oracle input whose original cost
is not included in this experiment; no deployable or end-to-end SAEM speed
claim follows. Loading and current preparation costs are reported.

## Methods and frozen budget

Compare the same baseline exact patient/population/shared-parameter updates
with and without recycling. Retain ordinary patient updates in both arms,
including 20% direct population proposals and pCN initialized with beta
sqrt(0.19). Use exact, ODE-free population-location updates and population-scale
updates, cached observation-noise updates, and exact dynamic-shared updates.
Adapt proposal scales only during the first 500 sweeps; retain draws afterwards.

One recycling pool has 36 slots for the nine ODE-requiring patients; equilibrium-
only patients continue their ordinary exact updates. Refresh four unused slots
per sweep. Assignments are injective. q is frozen; current population densities
are not. Common-grid proposal scores use the declared finite log floor -1e6,
with the exact MH correction against the original likelihood. Failed common
solves must not silently remove exact-positive support. Proposal cache changes
and exact accepted/rejected state changes must follow the tested kernel contract.

Both arms attempt the dynamic block every five sweeps, initially with log-scale
standard deviations (0.04, 0.04). Initial noise log-proposal SD is 0.06. Accepted
dynamic changes invalidate shared trajectory caches; rebuilding is counted.

Each chain has a maximum of 100,000 cost units and 20,000 sweeps. One unit is an
exact prediction attempt for an ODE-requiring patient, including a rejection
before integration, or a common-grid prediction attempt. Equilibrium-only
calculations consume zero such units but remain included in measured CPU/wall
time. Count actual integrations, failures, setup, warmup, reserve creation,
refresh, exact corrections, and dynamic-cache rebuilding separately. Common
solves are not assumed to cost the same CPU time as patient solves. No automatic
budget extension, second pool-size search, 115-patient chain, SAEM refit, neural
fit, or population-structure run follows this experiment.

Run four independent one-core workers concurrently, eight chains total. Save
checkpoints every 250 sweeps. An implementation smoke run uses one chain pair,
budget 1,500, warmup 10 and at most 60 sweeps; it cannot pass the scientific gate.

## Analysis and decision

Save states at completed sweeps on their regular iteration index, not at equal
cost timestamps. For each method use the longest equal-length retained prefix
across its four chains for joint diagnostics. Charge the full measured chain
cost, including unused tails, and report their lengths. Use the established R
posterior package rather than a new ESS/R-hat implementation.

For all 21 global coordinates and treated_logit_location = mu_u_eta_pi +
beta_nelf, require rank-normalized split R-hat <= 1.01, bulk ESS >= 400 and
tail ESS >= 200 in both arms. Compare posterior means, SDs and 5%, 50%, 95%
quantiles using the independent arms' combined Monte Carlo standard errors.
Require agreement within simultaneous 95% Bonferroni intervals over all these
comparisons. Report raw differences and standard errors; agreement alone does
not certify convergence or rule out a shared implementation defect.

The primary population summaries are beta_nelf, treated_logit_location,
log_omega_eta_pi, mu_log_lambda and log_omega_lambda. For bulk ESS per total
measured CPU time, require their median recycling/baseline ratio >= 1.5 and
minimum ratio >= 1.0; require each dynamic shared parameter's ratio >= 0.8.
Report the complete bulk/tail ESS-per-CPU and ESS-per-wall summaries, not only
favourable coordinates. Equal preparation cost is allocated to both arms;
the original comparator fit remains excluded and explicitly disclosed.

A pass supports this oracle-initialized reduced sampler comparison only. It
does not establish model-structure reuse, full-115 performance, or SAEM-level
end-to-end cost. Failure of convergence makes the accuracy/efficiency comparison
inconclusive; it must not be relabelled an impossibility result for recycling.
No failed gate triggers additional experiments automatically.
