# System A SAEM-message review packet

Date: 2026-09-04

This packet is intended for an independent review of one decision:

> Is it scientifically worthwhile to run one final bounded 12-patient test of
> the claim that a SAEM fit can supply sufficiently good patient anchors that
> only a short post-SAEM MCMC refresh is needed for useful fixed-psi Bayesian
> population uncertainty?

It is **not** evidence that the proposed method already works, nor is it a
self-contained reproduction archive: the exact aggregate evidence is included,
but the raw patient-state banks and fitted-object RDS files are not published.

## Current verdict

Do not scale to a 115-patient modular posterior. The last bounded pilot is now
complete and failed at both corrected SAEM branches. The status is no-go for
both the defensive-mixture implementation and the short-refresh pure-SAEM
single-anchor proposal.

The raw identity

```text
E_{X ~ L_i rho_i / Z_i}[g_i(X | eta, z_i) / rho_i(X)]
  = m_i(eta) / Z_i
```

is the existing Stage-1/Stage-2 message identity. The only potentially new
claim is operational: using SAEM for `rho_i`, starts, or branch anchors might
make patient-bank construction cheap enough to approach SAEM cost. That claim
has not yet been tested under a strict short-refresh budget.

## Corrected SAEM jobs

| Job | Seed | State | Elapsed |
|---|---:|---|---:|
| `122128` | 8 | completed, exit 0 | 12:11:34 |
| `122129` | 29 | completed, exit 0 | 12:01:03 |

The omega-initialization correction worked exactly in both jobs, but the final
fits are incompatible. In particular, implied untreated/treated PI efficacy
is approximately `0.99944 / 0.01919` for seed 8 and `1 / 0.999992` for seed
29. Neither fit is a defensible sole anchor, and no marginal likelihood was
computed to rank the modes. Do not select or average them after inspection.

The reviewable aggregate artifacts are under:

- `outputs/system_a_corrected_full_seed8_20260903a/`
- `outputs/system_a_corrected_full_seed29_20260903a/`

Raw fit RDS files, patient-level estimates, input observations, and the very
large repeated ODE-failure ledgers are intentionally not published.

## Defensive-reference calibration

Crucially, this calibration used the **old pre-correction production SAEM
anchor**. It did not test either of the two corrected seed-8 and seed-29 SAEM
branches reported above.

| Job | Role | State |
|---|---|---|
| `122130` | 12 defensive-reference patient banks | all 12 completed, exit 0 |
| `122142` | ODE-free pilot/plan | completed, exit 0 |
| `122144` | 24 component-target banks | all 24 completed, exit 0 |
| `122168` | initial bridge assessment | failed closed on nonconvergence |
| `122169` | structured-failure test suite | completed, exit 0 |
| `122170` | recorded bridge assessment | completed, exit 0 |

The tested reference was

```text
h_i = 0.5 g_SAEM,i + 0.5 g_priorcentral,i.
```

It is unusable as one patient MCMC bank for these data. Important diagnostics:

- mean sampling acceptance 9.39%, with 24/48 chains below 2%;
- minimum coordinate ESS 2.80;
- prior-component minimum relative weight ESS 0.001;
- prior-component pooled bridge nonconvergence for patients 3, 6, 55, and 74;
- prior-component cohort forward-chain range 1774.47;
- SAEM-component cohort forward/reverse log ratios 1.757 and 8.318.

Bank construction used 288,144 exact prediction calls and 213,709 successful
ODE integrations in total. The pure-SAEM component chains had much better mean
acceptance (42.175%) than the defensive reference, but their minimum coordinate
ESS was still only 6.44; that is evidence of unresolved local mixing, not a
successful pure-SAEM message test.

This is a no-go for the defensive mixture, not a test of nearby population
changes from a pure SAEM-conditioned bank. The exact CSV evidence is under:

- `outputs/system_a_message_plan_20260903a/`
- `outputs/system_a_message_assess_components_20260903a/`
- the `summary.csv` and `manifest.csv` files in the two bank directories.

## Final fork (completed)

The final single-anchor gate used the following predeclared design:

1. Treat both corrected SAEM endpoints as predeclared stress branches; do not
   choose the more convenient one.
2. Validate both against the sealed VODE-BDF target.
3. On the same 12 patients, target `L_i * g_SAEM,i` directly, without the
   defensive mixture.
4. Impose a short post-SAEM ODE budget in advance.
5. Reweight to genuinely nearby fixed-psi treatment-effect and population-scale
   endpoints.
6. Require patientwise MCMC diagnostics, importance-weight ESS/D2, independent
   candidate-bank/bridge agreement, and total cost accounting.

Jobs `122175` and `122176` completed the two pure-SAEM 12-patient banks, but
both ODE-free Stage-1 jobs (`122178`, `122177`) failed. Every endpoint failed.
Raw relative weight ESS remained above 0.20, but relevant pooled bulk ESS was
only 7.29--23.30 and 8.79--18.29, maximum split R-hat was 1.175 and 1.143, and
projected 115-patient treatment-message MCSE was approximately 0.90--1.01.
The candidate and bridge stage was therefore not run.

This closes the proposed fork under its own stop rule. The limiting issue is
not the algebraic identity or immediate local overlap; it is obtaining
reliable conditional patient draws at the claimed short post-SAEM cost. See
`FINAL_PURE_SAEM_RESULT.md` for the complete job and endpoint tables.

See `RESULTS.md` for the complete chronological interpretation and the source
files/tests for the exact estimators and failure semantics.
