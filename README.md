# PSAA: A physics informed neural network approach to quantify antigen presentation activities at single cell level using omics data
**PSAA** is an R package that offers a detailed, stepwise quantification of MHC pathway activity, enabling predictions of gene-specific impacts and their downstream effects on immune interactions. **PSAA** works with diverse omics data types, such as bulk data, single cell RNA-seq data, spatial transcriptomic data and proteomics data.

**PSAA** built upon [scFEA](https://github.com/changwn/scFEA) for robust metabolic flux estimation.

Briefly, **PSAA** takes as input: 1) a gene expression matrix (or protein expression matrix) stored in a Seurat object and 2) a MHC-I or MHC-II antigen presentation pathway, which has already been reconstructed in this paper, see details form []. The reconstructed pathways were represented as directed factor graphs where each intermediate state of antigen or MHC complexes is a factor and each reaction step is a variable. **PSAA** can further estimate the activity level of each reaction step using expression changes of the genes or proteins involved in the steps.
