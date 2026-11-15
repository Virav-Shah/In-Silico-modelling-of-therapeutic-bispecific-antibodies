# Working optimisation prototype

We now have a functioning inference-led prototype that finds high-tumor/low-healthy-kill designs in all three scenarios. It uses the same mechanistic workflow as the paper, with explicit inferred physiology and a broader engineering search space. It does not claim the original source-code implementation or a global optimum.

## Selected results at 24 hours

| Scenario | Tumor kill | Healthy kill |
|---|---:|---:|
| equal | 96.815% | 0.925% |
| unequal | 97.277% | 0.922% |
| low | 87.475% | 0.941% |

Healthy-cell kill is weighted over all healthy cells. Search feasibility requires <=0.95%, leaving 0.05 percentage-point nominal headroom below 1%. Selection maximises tumor kill among feasible evaluated candidates, after tighter-tolerance re-integration.

| Scenario | kon1 (M^-1 s^-1) | koff1 (s^-1) | kon2 (M^-1 s^-1) | koff2 (s^-1) |
|---|---:|---:|---:|---:|
| equal | 1640.14 | 1 | 982.378 | 0.899547 |
| unequal | 930.129 | 4.03154 | 100.47 | 0.615963 |
| low | 991.062 | 0.0969361 | 173.417 | 0.012768 |

## What we inferred and why

1. **Effective accessible cell burden**, rather than a guessed total-body cell count, controls the systemic binding sink. Fit value: approximately 4.73e12 reference cells.
2. **Effective intracellular volume** controls uptake-to-concentration conversion. Fit value: approximately 0.413 pL per cell, retaining the supplement's 500 µm² area. The paper never required a perfect spherical cell; real membrane area can exceed the area of a volume-equivalent sphere. This is an effective model parameter, not a measured anatomical value.
3. These two values were fitted only to approximate 24-hour tumor ADC concentrations read from Figure 4b (~56 µM) and Figure 6b (~280 µM). All other nuisance settings remained fixed. Matching those two anchors is calibration, not independent validation, and the fitted values are conditional on the remaining model assumptions.
4. Cell volume does not feed back into the current PK/binding equations. The ratio of the two uptake endpoints therefore constrains cell burden independently of the shared 1/volume concentration factor; their absolute scale then constrains volume. This gives a mechanistic reason for fitting these two parameters.
5. **Healthy uptake plotting:** the Figure 4 healthy curves appear to show the concentration averaged over target-positive healthy cells. Our predicted means are ~9.66 µM for the bispecific (all-healthy mean divided by 0.21) and ~16.90 µM for the monospecific (divided by 0.11), close to the figure's ~9 and ~17 µM. This is an inferred reporting convention. Toxicity still averages the separate Emax responses over ALL healthy cells; it is unchanged by the plotting convention.
6. **Search domain:** the figure heatmaps are not declared optimiser bounds. We explored kon=1e2–1e7 M^-1 s^-1 and koff=1e-5–1e2 s^-1. This permits weak initial capture with rapid second-arm rebinding. These are mathematical engineering candidates; experimental attainability of the selected kinetics is not established. Several selected on-rates lie below the paper's displayed sweep range.

## Baseline checks without further fitting

At the paper's fixed Figure 6 rates, the inferred model gives 98.59% tumor kill, 1.769% healthy kill, and 33.2% terminal occupancy, broadly consistent with its displayed behavior. The ~280 µM uptake is a calibration anchor, and its associated Emax kill is a direct consequence, not an independent validation point. The unoptimised bispecific gives 93.33% / 14.82% kill. The unoptimised monospecific gives 94.41% / 8.89% kill and 67.57 µM uptake, compared with roughly 74 µM in the figure. This recovers the expected unoptimised comparator ordering.

The unequal-density Figure 8 candidate remains a mismatch: predicted terminal T2 occupancy ~82.2% versus visually ~66%, uptake ~275 versus ~240 µM, and healthy kill ~4.59%. We did not fit this away. The new optimiser searches different binding rates under a strict healthy-kill constraint.

## Optimisation implementation

NSGA-II in Platypus; 48 candidates per generation; two independent seeds, 576 algorithm evaluations per scenario per seed; four log10-rate variables. Objectives are (-tumor kill, healthy kill); constraint is healthy fraction <=0.0095. Each evaluation integrates PK, binding, receptor turnover and ADC uptake for 24 hours. The selected best is retained from the entire evaluated archive, not just the last population. The selected states and full input parameters are saved in `results.json`.

- equal: seed 8407: 96.815% tumor / 0.925% healthy; seed 2026: 96.433% tumor / 0.871% healthy
- unequal: seed 8407: 97.232% tumor / 0.911% healthy; seed 2026: 97.277% tumor / 0.922% healthy
- low: seed 8407: 86.527% tumor / 0.816% healthy; seed 2026: 87.475% tumor / 0.941% healthy

The repeated searches give similar useful designs but are not a convergence or global-optimality proof. The unequal-density candidate is near the lower kon2 bound, so further bound sensitivity is warranted before treating it as a final molecular design.

## Verification and sensitivity

- Eight existing model tests pass.
- Independent BDF and Radau integrations agree within 1e-6 in endpoint kill fraction for all selected candidates (actual differences <1e-9).
- Drug mass balance closes within 1e-6 relative error (actual selected errors <1e-12).
- All selected candidates remain <=0.95% healthy kill on tighter and independent integration.
- Twelve frozen-design checks per scenario explore cell burden, volume, turnover, target density and the alternative surface off-rate convention. No reoptimisation is used in these checks.
- equal: largest tested healthy kill 1.125% under `turnover_double`.
- unequal: largest tested healthy kill 1.104% under `turnover_double`.
- low: largest tested healthy kill 1.525% under `turnover_double`.

Thus the nominal designs meet the constraint, but they are not guaranteed to stay below 1% under every assumption change. `calibration_reading_sensitivity.csv` additionally refits the two nuisance parameters at four illustrative combinations of 50–62 µM and 250–310 µM plot anchors; those intervals are reading-sensitivity choices, not statistical confidence intervals.

## Run

```sh
.venv/bin/python replication/prototype.py --nfe 576 --seeds 8407 2026
.venv/bin/python replication/validate_prototype.py
.venv/bin/python replication/calibration_sensitivity.py
.venv/bin/python replication/report_prototype.py
```

The earlier reconstruction and results remain available under `output/simulation/`; these new results are under `output/prototype/`. `configuration.json` records calibration and bounds; `*_evaluations.csv` records all screened candidates; `*_selected_trajectory.csv` records full occupancy and uptake readouts; `*_sensitivity.csv` records frozen-design perturbations. No user input is required to run this prototype. Measured target densities, trafficking and feasible binding-rate ranges would make subsequent designs more specific.
