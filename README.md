# code_Joint_model_SMART

This repo consists of the R codes to simulate SMART data with survival outcomes, fit the joint model described in "Joint Modeling and Multiple Comparisons with the Best of data from a SMART with Survival Outcomes". The codes were written by Yan-Cheng Chao and Qui Tran. Both authors contributed equally to the R codes.

In particular:
- simSMART.R: Function simSMART() simulates 2 stage SMART design II with time-to-response and overall survival outcomes from 4 regimens.
- simSMART_two_censoring.R: Functions simSMART_resp_dependent_censoring() and simSMART_covariate_dependent_censoring() simulate 2 stage SMART design II with time-to-response and overall survival outcomes from 4 regimens and incorporate response-dependent censoring and covariate-dependent censoring, respectively.
- calH.R: Function calH() to estimate the 2 baseline hazards (time-to-response and time-to-death) non-parametrically. 
- loglik.R: Function loglik() to calculate the log-likelihood for any set of beta coefficients.
- simulation_tab1.R: This contains the code to estimate the beta coefficients of the joint model.
- relative_efficiency.R: Contain procedure to calculate bootstrapped standard error (SE) for survival estimates and calculate Relative Efficiency between the proposed joint model and previous methods.
- multiple_comparison_with_the_best.R: Procedure to pick set of best regimens (regimens that has highest survival rates at a particular time of interest) by conducting multiple comparison with the best.
- data_analysis.R: Contain procedure to perform data analysis on CALGB dataset and generate the figure 3 in the manuscript.
- calgbdata.xls: The dataset from the CALGB trial.
