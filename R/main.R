#' Antigen presentation level prediction using single cell/bulk RNA-seq or spatial transcriptomic data
#' @param seurat_object Seurat object containing a query data set - antigen presentation level prediction will be applied to this object
#' @param pathway Select a pathway you want to used to predict the anitgen presentation levels. Choosing from "mhc1","mhc2" or "self_define"
#' @param protein Whether to used proteomics data
#' @param sample_name The name of your input data
#' @return seurat object with predicted antigen presentation levels saved in metadata
#' @export
psaa<-function(seurat_object,
               pathway,
               sample_name,
               protein=FALSE,
               seed=123
               )
{
  expression_matrix <- GetAssayData(seurat_object, slot = "data")
  current_path <- getwd()
  # Check if the folder 'scFEA' exists under the current path
  scFEA_path <- file.path(current_path, "scFEA")
  if (file.exists(scFEA_path) && file.info(scFEA_path)$isdir) {
    message("The 'scFEA' folder exists.")
  } else {
    message("The 'scFEA' folder does not exist. Downloading scFEA...")
    system("git clone https://github.com/changwn/scFEA")
  }
  # prepare data
  write.csv(expression_matrix,paste0(current_path,'/scFEA/data/',paste0(sample_name,'.csv')))
  # prepare pathway
  pathway_folder <- system.file("pathway", package = "PSAA")
  csv_files <- list.files(pathway_folder, pattern = "\\.csv$", full.names = TRUE)
  file.copy(from = csv_files, to = paste0(current_path,'/scFEA/data/'), overwrite = TRUE)

  filename=c()
  if(pathway=="mhc1"){
    system(paste0("python3 scFEA/src/scFEA.py --input_dir scFEA/data --res_dir scFEA/output --test_file ",sample_name,".csv --moduleGene_file mhc1_module_genes.csv --stoichiometry_matrix mhc1_cmMat.csv --cName_file mhc1_cName.csv"))
  }
  if(pathway=="mhc2"){
    system("cd scFEA")
    system(paste0("python3 scFEA/src/scFEA.py --input_dir scFEA/data --res_dir scFEA/output --test_file ",sample_name,".csv --moduleGene_file mhc2_module_genes.csv --stoichiometry_matrix mhc2_cmMat.csv --cName_file mhc2_cName.csv"))
  }
  output<-read.csv(paste0("scFEA/output/",filename,".csv"))
  M<-c()
  if(pathway=="mhc1"){
    M<-c('M_1','M_2','M_3','M_4','M_5','M_6','M_7','M_8')
  }
  if(pathway=="mhc2"){
    M<-c('M_1','M_2','M_3','M_4','M_5','M_6','M_7')
  }
  addfeature<-function(m){
    seurat_object<-AddMetaData(
      object = seurat_object,
      metadata = output[[m]],
      col.name = m
    )
    return(seurat_object)
  }
  for (i in 1:length(M)) {
    seurat_object<-addfeature(M[i])
  }
  ##antigen mean
  antigen_mean<-c()
  if(pathway=="mhc1"){
    antigen_mean<-(output[['M_1']]+output[['M_2']]+output[['M_3']]+output[['M_4']]+output[['M_5']]+output[['M_6']]+output[['M_7']]+output[['M_8']])/8
  }
  if(pathway=="mhc2"){
    antigen_mean<-(output[['M_1']]+output[['M_2']]+output[['M_3']]+output[['M_4']]+output[['M_5']]+output[['M_6']]+output[['M_7']])/7
  }
  seurat_object<-AddMetaData(
    object = seurat_object,
    metadata = antigen_mean,
    col.name = "ant_mean"
  )
  ##antigen presentated
  antigen_pre<-c()
  if(pathway=="mhc1"){
    antigen_mean<-(output[['M_5']]+output[['M_6']]+output[['M_7']])/3
  }
  if(pathway=="mhc2"){
    antigen_mean<-(output[['M_6']]+output[['M_7']])/2
  }
  seurat_object<-AddMetaData(
    object = seurat_object,
    metadata = antigen_pre,
    col.name = "ant_pre"
  )
  return(seurat_object)

}
#' Plot the predicted antigen presentation levels for spatial transcriptomics data
#' @param seurat_object the input seurat object
#' @param pathway Select a pathway you want to used to predict the anitgen presentation levels. Choosing from "mhc1","mhc2" or "self_define"
#' @param mode The type of antigen presentation levels you want to plot. Choosing from "all", "avg" ,"pre" ,M_1" ,"M_2"..."M_8"
#' @return NULL
#' @export
plot_levels<-function(seurat_object,
                      pathway,
                      mode){
  # Plot and save image(M)
  draw2<-function(g){
    s<-SpatialFeaturePlot(seurat_object, features = g)
    print(s)
    # plot(1)
    #dev.copy(jpeg,filename=paste0("plot", g, '.jpg'))
    png(paste0("plot/", g, '.png'))
    print(s)
    dev.off()
  }

  M<-c()
  if(pathway=="mhc1"){
    M<-c('M_1','M_2','M_3','M_4','M_5','M_6','M_7','M_8')
  }
  if(pathway=="mhc2"){
    M<-c('M_1','M_2','M_3','M_4','M_5','M_6','M_7')
  }
  if(mode=="all"){
    sapply(M, draw2)
  }
  if(mode=="avg"){
    draw2("ant_mean")
  }
  if(mode=="pre"){
    draw2("ant_pre")
  }
  else{
    draw2(mode)
  }

}


