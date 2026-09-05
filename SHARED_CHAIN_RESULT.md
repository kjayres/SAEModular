# Reduced exact recycling comparison: no-go at the tested budget

Job `122242` completed successfully on 2026-09-05 in 1h 23m 30s. All eight
chains finished: four baseline and four recycling, each with a 100,000 forward-
attempt budget on the same twelve patients. All 17 population and four shared
parameters moved. Both methods used the original informative priors and sealed
numerical likelihood. The common trajectory was a proposal only, with exact
original-likelihood correction. The [frozen contract](SHARED_POOL_CHAIN_CONTRACT.md)
was not changed after inspecting these results.

## Main findings

| Diagnostic | Baseline | Recycling |
|---|---:|---:|
| Largest global rank-normalised R-hat | 1.114 | 1.194 |
| Smallest global bulk ESS | 27.58 | 15.35 |
| Saved sweeps per chain, including 500 warmup | 9,258 | 5,236–5,398 |
| All required convergence diagnostics passed | No | No |

The reported primary median bulk ESS/CPU ratio (recycling / baseline) was
**0.446**, below the predeclared 1.5 gate. The minimum primary ratio was 0.127;
the worse of the two dynamic shared-coordinate ratios was 0.404. Both arms'
nonconvergence makes these efficiency estimates descriptive, not a definitive
stationary speed ratio. The independent-arm posterior-agreement gate also failed;
this does not establish sampler bias when neither arm has adequately converged.

Recycling produced only **1,900 changes in 191,349 assignment updates**, about
0.99%. The remaining updates selected the current state. No proposed assignment
was rejected by the final exact numerical correction. This points to candidate
usefulness/self-selection, rather than that correction, as the observed obstacle.
Dynamic-parameter changes, local-state replacements and reserve refreshes also
limit the lifetime of reusable predictions. Cheap cross-patient evaluation by
itself did not translate into useful enough movement in this retained experiment.

The actual spent budgets were 99,990 per baseline chain and 99,944–99,952 per
recycling chain. Preparation used 14.515 CPU seconds. Warmup, preparation and
unused retained tails are charged in the reported efficiency; historical oracle
comparator fitting is excluded and no deployable end-to-end cost claim is made.
The ledgers separate forward attempts from actual integrations and cached or
equilibrium-only calculations.

## Decision

**Do not promote this implementation to the main modular Bayes project, expand
to 115 patients, or claim that exact recycling works well for System A.** The
bounded practical test is negative; posterior accuracy at convergence remains
unresolved. Neither an impossibility theorem for recycling nor a proved bias in
the exact kernel follows from this result.

Do not automatically extend the chains or rebuild their proposal distribution.
A future revisit would need a specific, separately justified way to improve the
very low useful-assignment rate without hiding its construction cost. The
residual-component feasibility screen is an independent idea, not a repair of
this failed comparison.

[Compact evidence](results/shared_pool_chain) includes posterior diagnostics,
agreement checks, efficiency, costs, counters and likelihood-call ledgers. Raw
chains remain in `outputs/shared_pool_chain_122242/`. These outputs retain global
traces and final patient states, not a stationary fixed-global patient bank.
