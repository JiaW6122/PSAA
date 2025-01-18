#' Antigen presentation level prediction using single cell/bulk RNA-seq or spatial transcriptomic data
#'
#'
#'
#' @param seurat_object Seurat object containing a query data set - antigen presentation level prediction will be applied to this object
#' @param pathway Select a pathway you want to used to predict the anitgen presentation levels. Choosing from "mhc1","mhc2" or "self_define"
#' @param protein Whether to used proteomics data
#' @param sample_name The name of your input data
psaa<-function(seurat_object,
               pathway,
               sample_name,
               protein=FALSE,
               seed=123,
               )
{
  normalized_matrix <- GetAssayData(seurat_object, slot = "data")
  current_path <- getwd()
  write.csv(counts,paste0(current_path,'scFEA/data/',paste0(sample_name,'.csv')))
  if(pathway=="mhc1"){
    system("cd scFEA")
    system(paste0("python3 src/scFEA.py --input_dir data --res_dir output --test_file ",sample_name,".csv --moduleGene_file mhc1_module_genes.csv --stoichiometry_matrix mhc1_cmMat.csv --cName_file mhc1_cName.csv"))
  }
  if(pathway=="mhc2"){
    system("cd scFEA")
    system(paste0("python3 src/scFEA.py --input_dir data --res_dir output --test_file ",sample_name,".csv --moduleGene_file mhc2_module_genes.csv --stoichiometry_matrix mhc2_cmMat.csv --cName_file mhc2_cName.csv"))
  }
}
