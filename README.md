# TNBC_RF
This repository contains the code and resources used in the study titled "A 29-Gene Expression-Based Random Forest Model for Chemoresistance Prediction in Triple-Negative Breast Cancer". It includes scripts for preprocessing microarray gene expression data, performing differential expression analysis, and training classification models.

The primary datasets analyzed are GSE25066 and GSE20194. From these analyses, candidate genes were selected to be used as features for the Random Forest model.

It also includes:
- Code for normalizing the expression levels of all samples from the combined datasets, using .CEL files obtained from GEO.
- A file containing the normalized expression values for all samples (fRMA method).
- A script for differential expression analysis applied to the unified datasets.
- Text files listing the selected samples from each dataset, based on exclusion criteria.
- A file mapping probe (transcript) IDs to their corresponding gene symbols.
