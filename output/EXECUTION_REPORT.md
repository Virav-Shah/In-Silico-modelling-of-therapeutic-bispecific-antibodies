# Execution report

Status: workflow reconstruction; exact quantitative replication remains unresolved.

- Four published-rate scenarios simulated; all eight scientific consistency tests passed.
- Seven 5×5 parameter grids completed (175 evaluations), with efficacy, toxicity and ratio CSVs and plots.
- Two Platypus NSGA-II searches, seed 8407, 240 algorithm evaluations each; 234 and 237 evaluator calls respectively (reused solutions can avoid fresh calls).
- Initialisation: 16 log-domain corners, two specified seeds, random fill to population 40.
- Unequal-density search: no feasible returned candidate. Minimum observed healthy kill was 1.103013%. This is not a proof of infeasibility.
- Low-density search: best feasible returned tumor kill 55.522312%, healthy kill 0.916205%. Full selected inputs are in its status JSON. The selected candidate was re-integrated at tighter tolerance and remained feasible.
- 16-page LaTeX PDF compiled and visually reviewed; final compilation has no overfull-box or package warnings.

These results differ from the paper. Missing PK, cell geometry/burden and mechanistic conventions were not calibrated to force agreement. The project is ready for author-source parameter substitution and quantitative validation, not for a claim of exact reproduction.
