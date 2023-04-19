# Integrative Gene Regulatory Inference via two regression-based methods
---

This repository contains code and data gene-specific integration of binding sites to expression in regression-based Gene Regulatory Network inference in *Arabidopsis thaliana*.

In particular, the file `weightedRF_weightedLASSOtutorial.Rmd` is a notebook showing how to apply the weightedRF and weightedLASSO functions to the expression dataset of the root response to nitrate induction in *Arabidopsis thaliana*, a dataset from [Varala et al, 2018](https://www.pnas.org/doi/abs/10.1073/pnas.1721487115).

GRNs are inferred from expression profiles and binding sites information of the regulators into the promoter regions.

This notebook also shows how to compute precision and recall of the inferred GRNs on the DAPSeq data and how to compute prediction errors (MSE) for both models.


# Requirements
---

Before being able to use weightedRF, the modified C++ Random Forest implementation have to be set up on your computer. To do so, run in your system terminal from the root of this project :

```
R CMD SHLIB inference_functions/Cpp_dependencies/rfutils.c
R CMD SHLIB inference_functions/Cpp_dependencies/regTree.c
R CMD SHLIB inference_functions/Cpp_dependencies/regrf_mse.c
```

The R packages that you will need to install are :
doParallel, parallel, foreach, doRNG, randomForest, igraph, boot, stringr, tictoc, GENIE3
