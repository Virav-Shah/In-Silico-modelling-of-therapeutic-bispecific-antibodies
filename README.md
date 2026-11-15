# In-silico modelling of therapeutic bispecific antibodies

I built this project to investigate how antibody binding kinetics, target density and payload potency affect the selectivity of bispecific antibody–drug conjugates (bsADCs). It reconstructs the modelling workflow in [Stoddard et al. (2026)](https://doi.org/10.20935/AcadDrug8407), then searches for density-specific binding rates that maximise tumour kill while limiting healthy-cell kill.

## Implementation

The current computational notebook is **[Bispecific_ADC_Modelling_v2.ipynb](Bispecific_ADC_Modelling_v2.ipynb)**: definitions and assumptions, complete ODEs, model code, calibration, binding-rate and density sweeps, multiobjective optimisation, selected rates, results tables and plots. The notebook runs independently of saved results or external model scripts. Version 2 keeps the modelling theory, fixes equation rendering and uses local colour scales with contours for the binding landscapes. It contains no software-version table or report-generation code. The **[original notebook](Bispecific_ADC_Modelling.ipynb)** is preserved unchanged.

The model couples two-compartment pharmacokinetics to sequential target binding, receptor turnover, internalisation and an Emax payload response. Two effective physiological parameters are calibrated once to published uptake landmarks and then held fixed. NSGA-II explores the tumour/healthy-kill trade-off, followed by constrained local refinement, across nine density/potency scenarios. Each scenario selects its own four binding rates, with a nominal healthy-kill limit of 0.95% at 24 hours.

The binding framework draws on [Sengers et al. (2016)](https://doi.org/10.1080/19420862.2016.1178437), and compartment parameters on [Mould et al. (2007)](https://doi.org/10.1111/j.1365-2125.2007.02914.x). Additional modelling choices are explicitly identified in the notebook and report.

## Results and interpretation

Reoptimisation produces approximately 91–99% tumour kill at moderate/high densities under the nominal healthy-cell constraint. At the lowest density, the 4 µM payload achieves approximately 18%; changing to a 10 nM payload raises this to approximately 86%. Transferring a design between densities can violate the healthy-cell constraint, motivating separate optimisation.

This is a **reference-informed reconstruction, not an exact quantitative reproduction**. The two fitted uptake values are calibration targets. The unfitted Figure 8 comparison predicts roughly 277 µM uptake versus 240 µM in the paper, and 92% T2 occupancy versus 66%. Selected rates approach search bounds, and some physiological perturbations exceed 1% healthy kill. The results establish a working computational design workflow, not unique global optima, experimentally reachable antibodies or clinical safety.

See **[the PDF report](report/report.pdf)** for the equations, optimisation method, results and paper comparison. Its **[LaTeX source](report/report.tex)** and the figures it uses are included; the original notebook regenerates them.

## Run

Use Python 3.12 or newer in an environment with the dependencies in `requirements.txt` (tested with Python 3.14).

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter lab Bispecific_ADC_Modelling_v2.ipynb
```

Select **Restart Kernel and Run All Cells**. The complete run performs thousands of ODE integrations and takes several minutes or longer. Set `WORKERS = 1` in the first code cell for sequential execution. Search seeds and bounds are explicit; small numerical differences may occur across platforms.

For execution without the Jupyter interface:

```bash
jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=3600 Bispecific_ADC_Modelling_v2.ipynb
```

Run from the repository directory. Version 2 displays its results inline. To regenerate the existing report, run the original notebook with `pdflatex` installed.
