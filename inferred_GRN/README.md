# Final GRN optimized using DIOgene

We make available to the community the end product of DIOgene on the nitrate induction case study, for anyone interested in mining the inferred GRN.

This rdata file contains a list of two elements:

+ `"importances"` : a list of two matrices, that are the fully connected n_TFs*n_genes importance matrices for weightedLASSO and weightedRF (averaged accross 10 replicates)
+ `"grns"` : a list of 6 dataframes, that are the inferred edges from the importances of weightedLASSO and weightedRF and for three network densities (0.005,0.01,0.05).


The code to used to save these results can be found at the end of the `Gene_specific_optimisation_of_TFBM_integration.Rmd` notebook.

Don't hesitate to ask if you would like to access other intermediary data.