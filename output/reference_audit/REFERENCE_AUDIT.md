# Reference audit: recoverable parameters versus ADC-specific extensions

Checked 2026-11-15. This audit corrects the earlier blanket description of these inputs as missing. Values explicitly stated in an antecedent model should be treated as reference-based defaults. Whether an extension retained every value is a separate provenance question.

## Ref. 15: Sengers 2016

[Primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4968105/), equations 1–10, Tables 1–2.

- Three binding ODEs and algebraic receptor conservation are provided.
- Surface on-rates: 2.4e13 dm²/(mol s); current implementation matches.
- Surface off-rate mapping: q2 = koff2, q1 = koff1 * kon2 / kon1 for equal surface on-rates; current implementation matches.
- Spherical diameter 10 µm implies area 314.159 µm² and volume 0.523599 pL (geometrical derivation). Geometry is recoverable; whole-cell volume as payload distribution volume remains an extension.
- Constant solution antibody; internalisation excluded. Thus this reference alone cannot specify the complete ADC PK/internalisation model.

The linked supplementary ZIP returned a download challenge instead of an archive; it was not inspected. Published equations were accessible.

## Ref. 17: Mould 2007

[Primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2000651/), PK Results, Table 2.

Vc = 11.3 L, Vp = 41.5 L, Q = 1.05 L/h. These match current code. Elimination is Michaelis–Menten:

Jelim = Vmax(W) Cc / (Km + Cc)

Vmax = 1020 (W / 1e10 cells/L)^0.194 µg/h; Km = 338 µg/L.

The current illustrative CL = 0.1 L/h does not implement this reference. W is white-blood-cell concentration, not total-body cell number. Mould also supplies WBC dynamics; transplanting those dynamics would introduce a different PD model from the ADC paper's payload-kill readout.

## Refs. 18–19: Lee 2018 / 2019

[2018 primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5768681/) and [2019 author manuscript](https://escholarship.org/content/qt4w49b3fd/qt4w49b3fd_noSplash_61979d946bff17513d601b35fde73529.pdf).

These are experimental internalisation studies. They support dependence on antigen identity and expression ratio, rather than furnishing the ADC paper's entire ODE parameterisation. Their assay cell counts do not establish an in-vivo accessible cell burden.

## Ref. 20: Dua 2015

[TMDD tutorial](https://pmc.ncbi.nlm.nih.gov/articles/PMC4505827/).

The general framework supports target synthesis/degradation, complex internalisation and free-drug distribution. This provides a basis for the synthesis closure ksyn = kdeg R0 under baseline steady state, rather than making that closure an arbitrary invention. It does not identify the numerical cell burden for this study.

## Consequences for our implementation

1. Binding kinetics and PK distribution already have reference support; they should not be listed as wholly unavailable.
2. Clearance is a concrete discrepancy between our implementation and the referenced PK model. A reference-faithful comparison should use Michaelis–Menten elimination with an explicit WBC scenario.
3. Use Sengers geometry as a prespecified alternative to our fitted volume/500 µm² area. Do not claim the fitted pair is reference-derived.
4. Retain the ADC paper's stated densities, turnover, dose, potency and population proportions when combining the models; the reference experiments address different targets and contexts.
5. The coupled model must account for drug removal once: Mould's empirical target-related elimination and explicit target binding/internalisation may overlap. This is a modelling inference requiring an explicit combination rule, not evidence that either source is incorrect.
6. Neither the recovered binding model nor WBC concentration uniquely determines the effective numbers of tumour and healthy cells represented in the ADC simulation. Nor do they determine the ADC optimiser's bounds or selection rule.

No new simulation or claim of improved reproduction is made by this reference audit. Earlier results remain results of the documented reconstruction, not of the full Mould clearance model.
