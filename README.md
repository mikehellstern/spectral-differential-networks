# Spectral Differential Network Analysis for High-Dimensional Time Series

This repository contains code for "Spectral Differential Network Analysis for High-Dimensional Time Series" (AISTATS 2025).

## Abstract

Spectral networks derived from multivariate time series data arise in many domains, from brain science to Earth science. Often, it is of interest to study how these networks change under different conditions. For instance, to better understand epilepsy, it would be interesting to capture the changes in the brain connectivity network as a patient experiences a seizure, using electroencephalography data. A common approach relies on estimating the networks in each condition and calculating their difference. Such estimates may behave poorly in high dimensions as the networks themselves may not be sparse in structure while their difference may be. We build upon this observation to develop an estimator of the difference in inverse spectral densities across two conditions. Using an $\ell_1$ penalty on the difference, consistency is established by only requiring the difference to be sparse. We illustrate the method on synthetic data experiments and on experiments with electroencephalography data.


## Code

The structure of the code used in the paper is below. Note that the due to size restrictions, neither the data nor the results are included in this repo. However, the code to generate any necessary data and run the simulations are included as well as any links to external data used as is the case in the EEG analysis.

### simulations

This folder contains all the code to run and analyze the simulations in Section 4 of the paper

### srm_eeg64 

This folder contains all the code to run and analyze the EEG data in Section 5 of the paper. 
