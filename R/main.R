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
  write.csv(expression_matrix,paste0(current_path,'scFEA/data/',paste0(sample_name,'.csv')))
  filename=c()
  if(pathway=="mhc1"){
    system("cd scFEA")
    filename=system(paste0("python3 src/scFEA.py --input_dir data --res_dir output --test_file ",sample_name,".csv --moduleGene_file mhc1_module_genes.csv --stoichiometry_matrix mhc1_cmMat.csv --cName_file mhc1_cName.csv"),intern = TRUE)
  }
  if(pathway=="mhc2"){
    system("cd scFEA")
    filename=system(paste0("python3 src/scFEA.py --input_dir data --res_dir output --test_file ",sample_name,".csv --moduleGene_file mhc2_module_genes.csv --stoichiometry_matrix mhc2_cmMat.csv --cName_file mhc2_cName.csv"),intern = TRUE)
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


