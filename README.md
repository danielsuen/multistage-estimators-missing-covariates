## Multistage Estimators for Missing Covariates and Incomplete Outcomes

This repository contains simulation code for the paper:
> Daniel Suen and Yen-Chi Chen. "Multistage Estimators for Missing Covariates and Incomplete Outcomes". 2021.

The link to the arXiv manuscript is [here](https://arxiv.org/abs/2111.02367).

## File Guide

**binarytreat-simulations**

| File | Description |
| --- | ----------- |
| causal_sample_size.R | script for showing estimation procedure with different sample sizes for binary treatment |
| causal_sensitivity.R | script for running sensitivity analysis on causal simulated data |
| functions_causalsim.R | functions for generating binary treatment data and running simulations |
| plotting_causal_sample_size.R| plotting script for causal sample size simulation |
| plotting_causal_sensitivity.R | plotting script for causal sensitivity analysis simulation |

**coxmodel-simulations**

| File | Description |
| --- | ----------- |
| cox_sim_sample_size.R | script for showing estimation procedure with different sample sizes for Cox model |
| cox_sim_sensitivity.R  | script for running sensitivity analysis on Cox model simulated data |
| functions_coxsim.R | functions for generating Cox model data and running simulations |
| plotting_cox_sample_size.R | plotting script for Cox sample size simulation |
| plotting_cox_sensitivity.R  | plotting script for Cox sensitivity analysis simulation |
