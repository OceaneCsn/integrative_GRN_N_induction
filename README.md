# Integrative Gene Regulatory Inference via two regression-based methods


This repository contains code and data for the manuscript "Integration of binding sites to expression improves regression-based Gene Regulatory Network inference in *Arabidopsis thaliana*", Cassan et al, 2022.

In particular, the file `bRF_LASSO-D3S_tutorial.Rmd` is a notebook showing how to apply the bRF and LASSO-D3S functions to the expression dataset of the root response to nitrate induction in *Arabidopsis thaliana*, a dataset from [Varala et al, 2018](https://www.pnas.org/doi/abs/10.1073/pnas.1721487115).

GRNs are inferred from expression profiles and binding sites information of the regulators into the promoter regions.

This notebook also shows how to compute precision and recall of the inferred GRNs on the validation database [ConnecTF](https://connectf.org/), and how to compute prediction errors (MSE) for both models.