# mv_spectral_confounder
This repository contains code and plots for "A spectral confounder adjustment for spatial regression with multiple exposures and outcomes." The paper can be found on [arXiv](https://arxiv.org/abs/2506.09325).

## Overview 

The scientific data management plan is described in the Word document `Prim_SDMP.docx`. All relevant code, except for the health outcome data that is proprietary, is included in this repository.

To reproduce data analysis, run `data-analysis-code/data_analysis_simulated.R`. All files needed for this script are provided. Code for compile the results and making some plots are included in the same script. Plots from this setup are provided in the folder `data-analysis-code/simulated_plots` for verification. 

Alternatively, use `reproducibility_materials/data_analysis/data_analysis.Rmd` to check reproducibility.

## Organization

### data-analysis-code

The code to prepare the data is in `data_preparation.R`. For those that do not have access to Medicare data, we provide an alternative script `data_preparation_synthetic.R` that generates random health outcomes, so the code can be tested. Both files include code that downloads SVI data directly through the R package `findSVI`. 

The data provided include 
* `data/outputs/US_ZCTA_ACSvars_2015.csv`
* `data/inputs/zcta_state_region.csv`
* `data/inputs/nhgis0015_ds172_2010_zcta.csv`

The shape file should be downloaded as a zip file from [here](https://www2.census.gov/geo/tiger/TIGER2023/ZCTA520/) and then extracted into a folder.

The original health outcome data from Medicare, `health_data.csv` is not available to the public.

### Data analysis workflow

#### With access to Medicare data 

After the R objects are created from `data_preparation.R`, `final_analyses.R` can be used to run the data analysis. Necessary functions are in `functions_analysis.R`. Then `final_results_plots.R` and `final_results_plots_with_USM.R` can be used to make figures. Some additional R scripts---`cross_validation.R`, `model_selection.R`, `tensor_summary.R`, `sensitivity_check.R`---are also provided. 


#### Without access to Medicare data 

Data with simulated outcomes can be created from `data-analysis-code/data_preparation_simulated.R`, but all RDS files are also provided in the folder `data/simulated`. Use `data-analysis-code/data_analysis_simulated.R` to run data analysis on a smaller subset (North Carolina only), check results, and make plots. Necessary functions are in `data-analysis-code/functions_analysis.R`. 

### Simulation code

Two R scripts---`simulation.R` and `simulation_nonlinear.R`---in the folder `simulation-code` provide code to run simulations. In the same folder, `submit.csh` is the script to run jobs in HPC. Necessary functions are in `functions_sim.R`, and `sim_result.Rmd` provides code for make figures and summaries for simulation results. 

R functions for simulation are included in `simulation-code/functions_sim.R`. The code to run simulations is included in `simulation.R`. The two matrices used to generate data, `B1` and `B2`, are placed in `simulation/inputs`.

To test the function, below is an example:

```{r}
source("functions_sim.R")
Qmat <- makeQ(20)
EI <- eigen(Qmat)
Gamma <- EI$vectors
W <- EI$values
B <- bSpline(W,df=10,intercept=T)
load("simulation/inputs/B1.RData")
load("simulation/inputs/B2.RData")
beta <- matrix(c(3,-3,2,3,5,4,-2,4,5,7,2,-3,3,4,6,1,-4,2,3,4,5,1,3,4,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow = 10, ncol = 5, byrow = T)/10
set <- create_set7(n=400, nz=10, E=10, R=5, Q=Qmat, sig2X=c(3,3,3,3,3,3,3,3,3), sig2Z=rep(16,10), sig2theta=c(2,2,2,2,2), tau2r=c(4,4,4,4,4), lambda0=.9999, lambda1=.9, lambda2=.9, betas=beta, B1=B1, B2=B2, rho1=.9, rho2=.9, n1=20, bXZ=1, bw=1)
Ystar <- t(Gamma)%*%set$Y
Xstar <- t(Gamma)%*%set$X
fit <- MCMC_model4_8_v7(Ystar, Xstar, Qmat = Qmat, E = 10, n = 400, R = 5, Q = 5, L=10, K=5, B = B, iters = 5000, burn = 1000, full_result = T, W=W)
```

## Citation 

The arXiv paper can be cited as:

Prim, S.-N., Guan, Y., Yang, S., Rappold, A. G., Hill, K. L., Tsai, W.-L., Keeler, C. and Reich, B. J. (2025) A Spectral Confounder Adjustment for Spatial Regression with Multiple Exposures and Outcomes. *arXiv*:2506.09325.


