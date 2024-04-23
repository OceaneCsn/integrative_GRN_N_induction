# DIOgene : Data Integration Optimization for regression-based gene regulatory retwork inference


This repository contains code and data for the manuscript "Gene-specific optimization of data integration improves regression-based Gene Regulatory Network inference in _Arabidopsis thaliana_".
GRNs are inferred from expression profiles of the root response to nitrate induction in *Arabidopsis thaliana*, from [Varala et al, 2018](https://www.pnas.org/doi/abs/10.1073/pnas.1721487115) and from TF binding motifs information as a prior (TFBMs). 

**Our work introduces a novel methodology to optimize the intensity of data integration in a hypothesis-driven and gene-specific manner.**


# DIOgene: usage


The main R function allows the inference of integrative GRNs with an otpimized contribution of the prior data to gene expression.

For each target gene, it optimizes the contribution of prior data (here TFBMs) to gene expression by maximizing the improvement of MSE (accuracy in predicting the target genes expression)  as a function of prior integration over a simulated null hypothesis. Depending on the specified `model`, DIOgene will either rely on linear regressions (weightedLASSO), or non-linear regressions (weightedRF).

Here is how to use it:

```r
source('inference_functions/DIOgene.R')

DIOgene_results <- DIOgene(counts, genes, tfs, prior_matrix,
                    model = "non-linear",
                    nrep = 100,  nCores = 45,
                    ALPHAS = seq(0,1, by = 0.1))
```

**Arguments:**

+  `counts` Expression matrix (genes in rownames, conditions/samples in columns)
+  `genes` Vector of genes (must be present in the rownames of counts) to be used in GRN inference as target genes
+  `tfs` vector of genes (must be present in the rownames of counts) that are transcriptional regulators to be used as predictors in the regressions for GRN inference
+  `prior_matrix` Prior matrix Pi, a score (*e.g* 0, 0.5, 1) for each TF-target interaction.
Can contain NAs for pairs that do not have a prior value available. NAs will be turned into neutral priors of 1/2.
+  `model` "linear" or "non-linear": the type of regression that must be performed. "linear"" is for weightedLASSO, "non-linear" is for weightedRF.
+  `nrep` Number of repetitions of each regression model at each $\alpha$ value, to estimate the MSE mean and dispersion on the true and null datasets. Must be large enough (*e.g* 100 by default).
+  `nCores` Number of cores for multithreading. (Not supported on Windows)
+  `quiet` Removing prints of individual model estimation times? 
+  `ALPHAS` Set of $\alpha$ values to be explored. Default/recommended value: from 0 to 1 with a step of 0.1.
+  `S` Number of bootstrapped samples for robust regression models estimation. If NULL (default), it is automatically set as in the article depending on the chosen regression model (2000 for weightedRF or 500 for weightedLASSO).

**Output**

It returns a list containing:

+ `mats`: the list of importance matrices for each $\alpha$ value and each replicate,
on the true and shuffled datasets
+ `mse`: the normalized MSE values for each gene, $\alpha$ value and replicate on the true and shuffled datasets
+ `edges`: the list of inferred edges for each $\alpha$ value, replicate, and three densities
on the true and shuffled datasets
+ `alphas_opt`: the list of optimal alpha values per gene.
+ `diogene`: the importance matrix of the gene-specific GRN with optimal alpha. It can be used as is (fully connected GRN), or pruned to a sparse GRN using the functions `weightedLASSO_network` or `weightedRF_network`, depending on the chosen model.


A reproducible example on a handful of genes in our demonstration dataset is provided in the notebook `DIOgene_tutorial.Rmd`.
It further shows how the prune the inferred importance matrix, and how to plot the optimal value of $\alpha$ for a given gene.


# Repository overview

The other notebooks or scripts are available to further explore our code and reproduce the article's analyses:

+ Raw expression data and prepared rdata files used as inputs are contained in the `data` and `rdata` folders.

+ Gene expression and TF binding motifs input datasets are prepared in the folders `expression_dataset_preparation` and `TFBM_dataset_preparation`, respectively.

+ The folder `inference_functions` contains the source code of various functions inferring GRNs, optimizing data integration, and evaluating inferred GRNs.

+ The file `weightedRF_weightedLASSO_tutorial.Rmd` is a notebook showing how to apply the weightedRF and weightedLASSO functions to the root nitrate induction dataset. It demonstrates the basic used of our re-implemented linear and non linear integrative GRN algorithms. It also shows how to compute some quality metrics used in the article, namely prediction error on unseen samples (MSE), or precision and recall of the inferred GRNs against DAP-Seq data.


+ `Gene_specific_optimisation_of_TFBM_integration.Rmd` is a notebook launching the computation of numerous models over a range of $\alpha$ values, on true data but also in our synthetic null dataset (and thus requires significant computing resources). It then computes the optimal value of TFBM integration strength based on the proposed procedure for each gene, and shows several TFBM integration scenarios for different gene examples. It stores the results in rdata format.


+  `Inferred_GRNs_analysis.Rmd` loads the results generated by `Gene_specific_optimisation_of_TFBM_integration.Rmd` and studies the inferred GRNs. In particular, it displays the behavior of effective data integration, MSE, TFBM support as $\alpha$ increases, and also displays precision and recall analyses of inferred GRNs against DAP-Seq interactions. It can be used to reproduce the main figures of the paper.

+ `Functional_analyses.Rmd` further studies the list of genes for which TFBM are consensually integrated or not, and the modelling of nitrate signalling pathways in inferred GRNs. It can be used to reproduce the supp figures and supp tables of the article.


+ In the `inferred_GRN` folder, we make available to the community the end product of DIOgene on the nitrate induction case study, for anyone interested in mining the inferred GRN.

+ `Comparison_to_existing_methods.Rmd` contains code to reproduce our comparison of weightedRF and weightedLASSO to state of the art competitor algorithms for integrative regression-based GRN inference.


# Requirements
---

R > 4.0.1.


The list of packages to install in order to run DIOgene is : "parallel, tidyverse, igraph, pROC, PRROC, reshape2, stringr, foreach, doRNG, dplyr, ggh4x, org.At.tair.db, AnnotationDbi, randomForest, doParallel, glmnet, tictoc, ppcor, tidyr, ggplot2, ggpubr".

Before being able to use weightedRF, the modified C++ Random Forest implementation from iRafNet have to be set up on your computer (only once). To do so, run the following commands in your system terminal from the root of this project :

```
R CMD SHLIB inference_functions/Cpp_dependencies/rfutils.c
R CMD SHLIB inference_functions/Cpp_dependencies/regTree.c
R CMD SHLIB inference_functions/Cpp_dependencies/regrf_mse.c
```

These dependencies will then be automatically called when running the weightedRF function.

For the demonstration dataset preparation, the tool DIANE was used that should be installed as follows :

```
library(remotes)
install_github("OceaneCsn/DIANE")
```

