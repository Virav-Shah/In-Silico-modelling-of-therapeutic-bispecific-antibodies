# Information needed for exact replication

This is a checklist, not a message sent to anyone.

1. Original Python source, full dependency versions and execution commands for Figures 4–8 and S1–S5.
2. Complete ODEs, state definitions, initial conditions and units.
3. Vc, Vp, Q, clearance values and elimination law; how the Mould model was adapted without double-counting target-mediated clearance.
4. Total cell number, tumor/healthy counts, cell area and cell volume; treatment of Table 1's 100.03% total.
5. Exact interpretation of twice-hourly turnover, free-target synthesis, bound-target replenishment, dual-complex uptake and recycling.
6. 2D association values and mapping of surface dissociation constants to the four varied solution-binding rates, including detailed-balance corrections.
7. Monospecific bivalent equations and all arm/bond statistical factors.
8. Intracellular ADC/payload conversion; aggregation order for healthy kill and mean intracellular concentration.
9. Solver/integrator, tolerances, max steps, handling of infusion discontinuity and simulation endpoint.
10. Platypus algorithm and version, objective vector, constraints, variable transforms/bounds, initial population, seed(s), population size, evaluation budget and final selection rule.
11. Fixed parameters and resolutions for all heatmaps; raw numerical trajectory and heatmap arrays.
12. All S5 binding rates and whether its optimiser bounds differed from the main analysis.
13. Correction of Table 2's Figure 4 scenario labels and the Figure 8 narrative's affinity direction.
14. Explanation of Figure 8's plotted healthy-kill curve relative to the stated <1% constraint.
15. Correct references for Table S1 target densities and cell area.
