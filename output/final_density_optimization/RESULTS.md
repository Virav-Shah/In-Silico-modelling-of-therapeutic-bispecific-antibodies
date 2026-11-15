# Final density-specific optimisation

Our calibrated model; nominal predictions, not an exact paper reproduction.

| T1 density | T2 density | EC50 (µM) | Tumour kill (%) | Healthy kill (%) |
|---|---|---|---|---|
| 1e-13 | 1e-13 | 4.0 | 17.58 | 0.728 |
| 1e-12 | 1e-12 | 4.0 | 58.18 | 0.950 |
| 1e-11 | 1e-11 | 4.0 | 91.05 | 0.950 |
| 1e-10 | 1e-10 | 4.0 | 97.49 | 0.950 |
| 1e-09 | 1e-09 | 4.0 | 98.74 | 0.950 |
| 1e-09 | 1e-10 | 4.0 | 98.26 | 0.950 |
| 1e-09 | 5e-11 | 4.0 | 97.74 | 0.950 |
| 1e-09 | 1e-11 | 4.0 | 93.16 | 0.950 |
| 1e-13 | 1e-13 | 0.01 | 85.83 | 0.950 |

The same physiological calibration is frozen across all densities. Bounds: kon 1e2–1e7 M^-1 s^-1; koff 1e-5–1e2 s^-1. Healthy design limit 0.95%. Selected points are best verified feasible candidates, not certified global optima.

See summary.csv for exact rates, dissociation constants, transfer comparisons, seed agreement and sensitivity maxima. See configuration.json for every assumption. Some perturbations can exceed the nominal healthy limit; refer to sensitivity.csv.

[PDF](final_density_optimization_report.pdf) · [LaTeX](final_density_optimization_report.tex)
