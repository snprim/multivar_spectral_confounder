Reproducible code for “A Spectral Confounder Adjustment for 
Spatial Regression with Multiple Exposures and Outcomes”

.
├── data_analysis/
│   ├── adj.mat_nc.rds (adjacency matrix for North Carolina)
│   ├── adj.mat_south.rds (adjacency matrix for the south region of US)
│   ├── data_analysis.pdf
│   ├── data_analysis.Rmd
│   ├── data_nc.csv
│   ├── data_south.csv
│   ├── fit_svi_K10_L10_stB.rds (saved results from MSM model with simulated outcomes for NC)
│   ├── fit_svi_naive.rds (saved results from naive model with simulated outcomes for NC)
│   └── functions_analysis.R
├── simulation/
│   ├── functions_sim.R
│   ├── sim_result.pdf
│   ├── sim_result.Rmd
│   ├── simulation.R (simulation settings for linear confounding and sensitivity analysis)
│   ├── simulation_nonlinear.R (simulation settings for nonlinear confounding)
│   └── submit.csh (a script to send the simulation code to HPC)
├── acc_form.pdf
├── acc_form.Rmd
└── Readme.txt


Workflow for rerunning data analysis:  All code is in data_analysis.Rmd. A saved output data_analysis.pdf is provided.

Relevant simulation code is provided. Simulation results can be found in sim_result.pdf. 