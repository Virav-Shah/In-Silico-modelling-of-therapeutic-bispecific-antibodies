# Implementation notes

The new density model deliberately selects `surface_off_rule='symmetric'`. Previous model profiles retain their original choice by default.

The unequal membrane-association-rate generalisation of the original Sengers closure contained an inverted ratio. With this code's labels, detailed balance requires:

q1 = koff1 * kon2 * k2d1 / (kon1 * k2d2), q2 = koff2.

The previous implementation interchanged k2d1 and k2d2 in this expression. Equal membrane on-rates cancel the ratio, so all previous Figure 4/6/8 and density-study endpoints are unchanged. The historical S3 diffusion sweep, which varies the two membrane rates independently, and its old generalised equation in model_analysis.pdf are superseded on this point. A corrected reference-physiology S3 sweep is included in corrected_S3/; its filename identifies a numerical correction, not new calibration or paper validation. The symmetric closure and the corrected general cycle are stated in the new report. Tests now use unequal membrane rates as well as unequal solution rates.

Numerical optimisation allows a 1e-7 kill-fraction tolerance (1e-5 percentage points) around the nominal 0.0095 boundary. Every selected candidate is independently verified below the strict paper threshold of 0.01.

The selected designs are nominal and bound-dependent. Some local searches hit their iteration limit; their best feasible evaluations are retained and the termination messages are saved. Independent seeds agree closely on achieved efficacy, but this does not prove a global or unique optimum.
