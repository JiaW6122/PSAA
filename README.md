# PSAA: A physics informed neural network approach to quantify antigen presentation activities at single cell level using omics data
**PSAA** is an R package that offers a detailed, stepwise quantification of MHC pathway activity, enabling predictions of gene-specific impacts and their downstream effects on immune interactions. **PSAA** works with diverse omics data types, such as bulk data, single cell RNA-seq data, spatial transcriptomic data and proteomics data.

**PSAA** built upon [scFEA](https://github.com/changwn/scFEA) and [scFEA-mpo](https://github.com/ptdang1001/scFEA/tree/main) for robust metabolic flux estimation.

Briefly, **PSAA** takes as input: 1) a gene expression matrix (or protein expression matrix) stored in a Seurat object and 2) a MHC-I or MHC-II antigen presentation pathway, which has already been reconstructed in this paper, see details from [Pathway Construction](https://github.com/JiaW6122/PSAA/blob/main/supplementary%20files/Pathway_Construction.md). The reconstructed pathways were represented as directed factor graphs where each intermediate state of antigen or MHC complexes is a factor and each reaction step is a variable. **PSAA** can further estimate the activity level of each reaction step using expression changes of the genes or proteins involved in the steps.

# Preparation
Install **PSAA** from GitHub:
```R
library(remotes)
remotes::install_github("JiaW6122/PSAA")
```

Set up the environment
```R
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reticulate)
``` 
```R
# Create a Conda environment
conda_create("scFEA_env")

# Install PyTorch
conda_install("scFEA_env", packages = c("pytorch", "torchvision", "torchaudio"), channel = "pytorch")

# Install magic-impute
conda_install("scFEA_env", packages = "pip")
py_install("magic-impute", envname = "scFEA_env")
conda_install("scFEA_env", packages = c("numpy",  "matplotlib"))
conda_install(envname = "scFEA_env", packages = "pandas")
```

```R
use_condaenv("scFEA_env", required = TRUE)
py_config()  # Verify the environment and Python path
```

# Tutorials 

## Spatial transcriptomics data


Download the human breast cancer data from: https://drive.google.com/file/d/1032ck3O4G72SnVSlEuuxt-m4vJIsOKgk/view?usp=share_link

Load demo dataset:
```R
HumanbreastV1 <- readRDS("your_path/HumanbreastV1.rds")
# HumanbreastV1 <- SCTransform(HumanbreastV1, assay = "Spatial", verbose = FALSE)
```

Run PSAA to predict sample-wise antigen presentation levels through MHC calss I antigen presentation pathway:
```R
library(PSAA)
HumanbreastV1 <- PSAA::psaa(HumanbreastV1, "mhc1", "Humanbreast")
```

Visualize the predicted antigen presentation levels in spatial transcriptomics data:
```R
PSAA::plot_levels(HumanbreastV1, "mhc1", "pre")
```

