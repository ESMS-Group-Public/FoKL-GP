# FoKL-GP
Karhunen Loève decomposed Gaussian processes with forward variable selection

FoKL-GP is a Gibbs sampler and automatic function builder for Karhunen-Loève decomposed Gaussian processes.
It uses a forward variable selection method as a way to contend with the dimensionality curse.
The routine works best for tabular data with a moderate number of continuous inputs, all of which should be 'min-max' scaled to the [0,1] interval.

Clearing outliers in the input dataset is also quite important.

The main routine is 'emulator' -- notes in the preamble to that function describe the inputs and outputs.

The coefficients ('betas') and interaction matrix ('mtx') -- which specifies which terms appear in the model -- may be used directly within the main GP evaluation routine 'bss_eval'.

The multiple-input, single-output GP functions are constructed through an automated model building routine utilizing linear Bayesian inference (through Gibbs sampling) on a sequence of models of increasing complexity. The Bayesian information criterion (or Akaike) is used to track model fitness scores. The routine requires a number of increases in the BIC (the BIC optimum is a minimum) before model termination, at which point the optimum model is returned.

There is a test case in the MATLAB workspace in the folder. Load the workspace and run 'emulator' with

[betas, mtx, evs] = emulator([reshape(X, m*n,1) reshape(Y, m*n, 1)], reshape(DATA_nois, m*n, 1), 0.01);

visualize it with

[meen, bounds, rmse] = coverage(betas, [reshape(X,m*n,1) reshape(Y,m*n,1)], [1 1], [], reshape(DATA_nois,m*n,1), mtx, 'standard', 100, 1);

Change the hyperparameters to see how various settings affect the routine.

A basis set from the BSS-ANOVA kernel is included -- if you use it please cite

Reich BJ, Storlie CB, Bondell HD. Variable selection in Bayesian smoothing spline ANOVA models: Application to deterministic computer codes. Technometrics. 2009 May 1;51(2):110-120. doi: 10.1198/TECH.2009.0013.

If you use FoKL-GP with any basis, please cite

Hayes, K, Fouts MW, Baheri A, Mebane DS. Fast Dynamic System Identification with Karhunen-Loève Decomposed Gaussian Processes. arXiv:2205.13676 

HOW TO USE OTHER GP KERNELS WITH FoKL-GP:

The BSS-ANOVA basis functions are built with spline_coefficient_500.txt and 'splineconvert500'.
Read the coefficients in the file into a single matrix and pass that to the splineconvert routine.
Other kernels can be used via the following procedure to generate the KL basis.

1. Use the kernel to generate a covariance matrix on a normalized [0,1]x[0,1] domain with 499 intervals.

2. Eigendecompose, ordering the normalized eigenfunctions by the absolute value of their eigenvalues. 

3. Scale each eigenfunction by the square root of its eigenvalue.

4. Fit to cubic splines where the input is normalized to each sub-interval.

5. Store the coefficients in a text file where the rows correspond to intervals and columns correspond to coefficients in ascending order. Stack the coefficients in the file without gaps, starting with the lowest order function.

6. Read the data into the MATLAB workspace and run 'splineconvert500'.


Work in this repository was supported by the National Science Foundation under an EPSCoR Track-2 Award (#2119688).
