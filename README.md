# Gene-specific integrative Gene Regulatory Network Inference for regression-based methods
---

This repository contains code and data for the manuscript "Gene-specific optimization of data integration improves regression-based Gene Regulatory Network inference in _Arabidopsis thaliana_".
GRNs are inferred from expression profiles of the root response to nitrate induction in *Arabidopsis thaliana*, from [Varala et al, 2018](https://www.pnas.org/doi/abs/10.1073/pnas.1721487115) and from TF binding motifs information (TFBMs). 

**Our work introduces a novel methodology to optimize the intensity of data integration in a hypothesis-driven and gene-specific manner.**

# Repository overview

+ Raw expression data and prepared rdata files used as inputs are contained in the `data` and `rdata` folders.

+ Gene expression and TF binding motifs input datasets are prepared in the folders `expression_dataset_preparation` and `TFBM_dataset_preparation`, respectively.

+ The folder `inference_functions` contains the source code of various functions inferring GRNs, optimizing data integration, and evaluating inferred GRNs.

+ The file `weightedRF_weightedLASSO_tutorial.Rmd` is a notebook showing how to apply the weightedRF and weightedLASSO functions to the root nitrate induction dataset. It demonstrates the basic used of our re-implemented linear and non linear integrative GRN algorithms. It also shows how to compute some quality metrics used in the article, namely prediction error on unseen samples (MSE), or precision and recall of the inferred GRNs against DAP-Seq data.


+ `weightedRF_weightedLASSO_gene_specific.Rmd` is a larger notebook launching the computation of numerous models over a range of $\alpha$ values, on true data but also in our synthetic null dataset (and thus requires significant computing resources). It then computes the optimal value of TFBM integration strength based on the proposed procedure for each gene and displays the results. It can be used to reproduce figures, supp figures and supp tables of the paper.


# Requirements
---

R > 4.0.1.


Before being able to use weightedRF, the modified C++ Random Forest implementation from iRafNet have to be set up on your computer. To do so, run the following commands in your system terminal from the root of this project :

```
R CMD SHLIB inference_functions/Cpp_dependencies/rfutils.c
R CMD SHLIB inference_functions/Cpp_dependencies/regTree.c
R CMD SHLIB inference_functions/Cpp_dependencies/regrf_mse.c
```

These dependencies will then be automatically called when running the weightedRF function.

Many packages are used in the R scripts, and can be installed from the CRAN in the usual way with install.packages(), except for DIANE, that should be installed as follows :

```
library(remotes)
install_github("OceaneCsn/DIANE")
```

