# Define QC thresholds
min_genes <- 200    # Minimum number of features per cell
max_genes <- 2500   # Maximum number of features per cell
max_mito <- 5       # Maximum percentage of mitochondrial genes

# Apply QC and normalization
ifnb.list <- lapply(ifnb.list, function(seurat_obj) {
  # Calculate percent.mt (percent mitochondrial genes)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # Subset cells based on QC thresholds
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mito)

  # Normalize the data
  seurat_obj <- NormalizeData(seurat_obj)

  # Return the processed object
  return(seurat_obj)
})

# Check the output
ifnb.list[[1]]
