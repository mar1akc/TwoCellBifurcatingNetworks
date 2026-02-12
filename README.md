# TwoCellBifurcatingNetworks
Codes for generating figures depicting the dynamics and bifurcations of two-cell feedforward networks of pitchfork and Stuart-Landau cells

### Folder PitchforkCells

**system1.m**

Set the variable flag in line 11:

* flag = 0: plot the phase diagrams in Figs. 5 and A1;
* flag = 1: plot the bifurcation diagrams in Fig. 5;
* flag = 2: plot the basins in Fig. 6;
* flag = 3: plot jumps in Fig. 6.

**System1_singularity.m** generates Fig. 8.

### Folder StuartLandauCells 

**system2_reduced.m** generates Figs. 10 and 11. To generate Fig. 11, uncomment line 215.

**system2_reduced_simulations.m** generates Fig. 12.

**System2_invariant_tori.m** generates Fig. 13. Set parameters in line 22.

**system3_reduced.m** generates Fig. 15.

**system3_reduced_mu_epsilon_paths.m** generates Fig. 16.

**system3_reduced_XPPAUT_plots.m** generates Fig. 17. To plot the insets, uncomment an appropriate axis command in lines 234—236.

**system4_reduced_gamma1.m** generates Fig. 19.

**system4_singularity.m** generates Fig. 22.

**system2_HBpts.m** generates Fig. 23 (left).
**system4_HBpts.m** generates Fig. 23 (right).

**system4_amplitude_comparison.m** generates Fig. 24.

**check_formulas.m** checks Eqs. (62) and (68).
