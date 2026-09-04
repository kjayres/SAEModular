# Shared-trajectory feasibility screen

Frozen before numerical results, 2026-09-04. This is a separate experiment from
the stopped Fable/SAEM message work. No retained sampler, surrogate training,
patient-bank construction, or SAEM rerun is authorized by this screen.

## Questions and design

1. Can actual reserve draws plausibly serve several System A patients?
2. Can one common-trajectory evaluation serve their observation designs
   accurately and materially faster than separate sealed evaluations?

Use the existing 12 audited patients, in canonical order:
3, 6, 16, 20, 55, 71, 74, 88, 105, 111, 117, 122.
Use two coherent Bayesian comparator snapshots: chains 1 and 3, retained row
750, with all population, shared and patient coordinates from that same row.
Artifact hashes and the mapping from patient IDs to Stan indices are checked.
These are oracle, posterior-relevant test contexts, not new conditional banks,
proof of comparator convergence, or deployable SAEM initializations.

Freeze one normalized reference q: the equally weighted mixture of the treated
and untreated population distributions at both snapshots (four components).
The use of Gaussian components reflects the present System A model only; the
pool construction and diagnostics accept arbitrary evaluable densities.
Generate 128 independent reserve vectors once, seed 20260904, and reuse these
same vectors under both shared-parameter snapshots. Do not redraw failed
equilibria or integrations. All remain in cost and usefulness denominators.

For each snapshot, evaluate these reserves and its 12 current patient vectors
against all 12 patients using both numerical interfaces. The current vectors
provide paired current likelihoods and a secondary cross-patient diagnostic.
They are occupied states, not draws from q; their cross-patient scores describe
potential future reuse after vacating a slot, not immediate legal assignments.
Exclude self-patient comparisons from that secondary summary.

Also audit the full 115-patient observation-grid union on eight preselected
vectors per snapshot: the first four reserves and the current states of panel
patients 3, 20, 74 and 122. This is a numerical-interface audit, NOT a 115-patient
inference run. Maximum prediction-interface calls: 5,496, at most one actual
integration per call. Count attempts, pre-integration rejections and integration
failures separately. No additional scientific solves are used to fit anything.

## Exact and experimental quantities

The reference is the unmodified, hash-pinned VODE-BDF target. At each snapshot:

    log w_i(z) = log L_i(z; psi) + log g_i(z; eta) - log q(z)
    alpha_i(z) = min(1, exp(log w_i(z) - log w_i(x_i)))

These are one-reserve MH probabilities. Report population-scale squared patient
jumps weighted by alpha, not population ESS or achieved chain mixing. The sum
of acceptances across patients is potential alternative usefulness: injective
assignments forbid multiple patients occupying one reserve simultaneously.

The common solve uses the sorted union of each patient's ALREADY adjusted
positive times, the same ODE/equilibrium/tolerances, and the original observation
likelihood. It has its own experimental numerical-target fingerprint. It is
not declared identical to, nor substituted into, the sealed numerical target.
Compare finite support, failure outcomes, log likelihoods and resulting MH
probabilities; no interpolation or solver-failure filtering is allowed.
Alternate which interface is timed first. Include likelihood reconstruction in
both timings, and report total and successful-integration-only timings. Use four
workers for independent vector evaluations; do not simulate parallel assignments.

## Predeclared operational decision

Both snapshots must satisfy all the following to justify discussing one reduced
sampler comparison (not to claim the method works):

- At least 10% of direct q reserves have alpha >= 0.1 for at least two patients.
- At least 6 of 12 patients have mean direct-reserve alpha >= 0.02.
- All paired cases agree on finite support and prediction success; maximum
  absolute finite log-likelihood difference <= 0.01, including the full-grid audit.
- Maximum absolute change in panel reserve MH probability <= 0.01.
- Separate/common aggregate elapsed-time ratio >= 3 for both the 12-patient
  panel and the 115-grid audit, with at least four common integrations per audit.

Report binomial Wilson intervals for the reserve-sharing fraction and all raw
patientwise metrics. These thresholds are practical screening choices, not
theorems or universal requirements. A failed gate stops THIS fixed-q/shared-grid
route under this design; it cannot rule out every possible q, solver or modular
Bayesian method. A pass still leaves dynamic-psi cache rebuilding, exclusion,
population coupling, preparation amortization and posterior agreement untested.

Save compact diagnostic CSVs, source/artifact hashes, complete call ledger and
an explicit result memorandum. Preserve older negative evidence. Do not delete
anything or start a neural/surrogate repair automatically.
