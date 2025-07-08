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


  # run scfea
  setwd(paste0(current_path,"/scFEA"))
  filename=c()
  if(pathway=="mhc1"){
    # system(paste0("python3 src/scFEA.py --input_dir data --res_dir output --test_file ",sample_name,".csv --moduleGene_file mhc1_module_genes.csv --stoichiometry_matrix mhc1_cmMat.csv --cName_file mhc1_cName.csv"))
    system(paste0("python src/scFEA.py --input_dir data --res_dir output --test_file ",sample_name,".csv --moduleGene_file mhc1_module_genes.csv --stoichiometry_matrix mhc1_cmMat.csv --cName_file mhc1_cName.csv"))
  }
  if(pathway=="mhc2"){
    system("cd scFEA")
    system(paste0("python src/scFEA.py --input_dir data --res_dir output --test_file ",sample_name,".csv --moduleGene_file mhc2_module_genes.csv --stoichiometry_matrix mhc2_cmMat.csv --cName_file mhc2_cName.csv"))
  }

  # read the result
  # Define the folder path
  output_folder <- "output"  # Replace with the actual path to your output folder

  # List all .csv files starting with sample_name
  csv_files <- list.files(
    path = output_folder,
    pattern = paste0("^",sample_name,".*\\.csv$"),
    full.names = TRUE
  )

  # Extract the timestamp from each file name
  timestamps <- sapply(basename(csv_files), function(file) {
    # Extract the timestamp using a regex
    match <- regmatches(file, regexpr("[0-9]{8}-[0-9]{6}", file))
    if (length(match) > 0) {
      return(match)  # Return the timestamp if found
    } else {
      return(NA)  # Return NA if no timestamp is found
    }
  })

  # Filter out files without timestamps
  valid_files <- !is.na(timestamps)
  csv_files <- csv_files[valid_files]
  timestamps <- timestamps[valid_files]

  # Convert timestamps to a datetime object for sorting
  timestamps <- as.POSIXct(timestamps, format = "%Y%m%d-%H%M%S")

  # Find the index of the newest file
  newest_index <- which.max(timestamps)

  # Get the newest file
  newest_file <- csv_files[newest_index]


  output<-read.csv(newest_file)
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
    antigen_pre<-(output[['M_5']]+output[['M_6']]+output[['M_7']])/3
  }
  if(pathway=="mhc2"){
    antigen_pre<-(output[['M_6']]+output[['M_7']])/2
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
#' @param sample_name The name of your input data
#' @param mode The type of antigen presentation levels you want to plot. Choosing from "all", "avg" ,"pre" ,M_1" ,"M_2"..."M_8"
#' @return NULL
#' @export
plot_levels<-function(seurat_object,
                      pathway,
                      sample_name,
                      mode){
  # Plot and save image(M)
  draw2<-function(g){
    s<-SpatialFeaturePlot(seurat_object, features = g)
    print(s)
    # plot(1)
    #dev.copy(jpeg,filename=paste0("plot", g, '.jpg'))
    png(paste0("plot/", sample_name,"_", g, '.png'))
    print(s)
    dev.off()
  }

  # Check if the "plot" folder exists
  setwd("..")
  if (!dir.exists("plot")) {
    # Create the folder if it doesn't exist
    dir.create("plot")
    message("The 'plot' folder has been created.")
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
    break
  }
  if(mode=="avg"){
    draw2("ant_mean")
    break
  }
  if(mode=="pre"){
    draw2("ant_pre")
    break
  }
  else{
    draw2(mode)
  }

}

#' Spatial segmentation upon varied dependency between predicted antigen presentation and T cell infiltration
#' @param seurat_object description
#' @param pathway description
#' @param sample_name description
#' @param k the number of neighbors we want to use
#' @param alpha the threshold of weight
#' @param sigma parameter of Gaussian kernel weight matrix
#' @return NULL
#' @export
spatial_seg<-function(seurat_object,
                      pathway,
                      sample_name,
                      k=100,
                      alpha=0.85,
                      sigma=10
                      ){
  # get gene expression using gene name
  tcell_gene<-c()
  if(pathway=="mhc1"){
    tcell_gene<-c("CCL5", "CD3D", "CD7",  "CD8A", "CST7", "CTSW", "GZMA", "GZMB", "GZMH", "GZMM", "NKG7", "PRF1")
  }
  if(pathway=="mhc2"){
    tcell_gene<-c("AIM2",     "ANKRD44",  "ARHGAP25", "ARHGAP9",  "ARHGEF6",  "BTLA",     "CD28",     "CD40LG",   "CD48",     "CD5",
                   "DOCK10",   "GIMAP1",   "GIMAP5",   "GPR171",   "GPR174",   "IKZF1" ,   "IRF8",     "ITK" ,     "LTA"  ,    "LY9"  ,
                   "PARP15",   "PTPRC",    "RASSF5",   "RHOH",     "SAMD3",    "SCML4",    "SLAMF1",   "TAGAP", "TESPA1" ,  "THEMIS",
                   "TMC8",     "TRAF3IP3", "TRAT1",    "ZC3H12D",  "ZNF831")
  }
  #compute T cell level of central point
  tcell_cen<-function(visium){
    # FetchData can pull anything from expression matrices, cell embeddings, or metadata
    tcell_level<-sapply(tcell_gene,function(x){return(FetchData(object = visium, vars =x))})
    tl<-matrix(unlist(tcell_level),ncol=length(tcell_gene))
    mean_level<-apply(tl, 1, mean)
    return(mean_level)
  }
  # find k nearest neighbors for each sample
  findk<-function(coords,k){
    dis<-distance(coords,method = 'euclidean')
    dd<-t(apply(dis,1,order))
    knestest<-dd[,2:k+1]
    return(knestest)
  }
  # compute Gaussian kernel weight matrix
  gua_kernel_weight<-function(coords,alpha){
    dis<-distance(coords,method = 'euclidean')
    kg<-exp(-dis/(2*sigma^2))
    # normailze
    # kg<-apply(kg,1,function(x){
    # x=x/sum(x)
    # return(x)
    # })
    kg2<-apply(kg,1,function(x){
      x[x<alpha]=0
      return(x)
    })
    # dd<-t(apply(dis,1,order))
    # knestest<-dd[,2:k+1]
    return(kg2)
  }
  tcell_level<-tcell_cen(seurat_object)
  # seurat_object<-AddMetaData(
  #   object = seurat_object,
  #   metadata = tcell_level,
  #   col.name = "tcell"
  # )
  coords = GetTissueCoordinates(seurat_object)
  weight<-gua_kernel_weight(coords,alpha)
  rownames(weight)<-NULL
  colnames(weight)<-NULL
  tcell_amp<-weight%*%array(tcell_level,dim=c(length(tcell_level),1))

  seurat_object<-AddMetaData(
    object = seurat_object,
    metadata = tcell_amp[,1],
    col.name = "tcell_amp"
  )
  antigen_pre<-FetchData(object = seurat_object, vars ="ant_pre")
  antigen_amp<-weight%*%array(antigen_pre,dim=c(length(tcell_level),1))
  seurat_object<-AddMetaData(
    object = seurat_object,
    metadata = antigen_amp[,1],
    col.name = "ant_amp"
  )

  #compute bivariate spatial correlation
  x<-antigen_amp
  y<-tcell_amp
  # Bivariate Moran's I
  moran_I <- function(x, y = NULL, W){
    if(is.null(y)) y = x

    xp <- scale(x)[, 1]
    yp <- scale(y)[, 1]
    W[which(is.na(W))] <- 0
    n <- nrow(W)

    global <- (xp%*%W%*%yp)/(n - 1)
    local  <- (xp*W%*%yp)

    list(global = global, local  = as.numeric(local))
  }
  # Permutations for the Bivariate Moran's I
  simula_moran <- function(x, y = NULL, W, nsims = 10000){

    if(is.null(y)) y = x

    n   = nrow(W)
    IDs = 1:n

    xp <- scale(x)[, 1]
    W[which(is.na(W))] <- 0

    global_sims = NULL
    local_sims  = matrix(NA, nrow = n, ncol=nsims)

    ID_sample = sample(IDs, size = n*nsims, replace = TRUE)

    y_s = y[ID_sample]
    y_s = matrix(y_s, nrow = n, ncol = nsims)
    y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)

    global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
    local_sims  <- (xp*W%*%y_s)

    list(global_sims = global_sims,
         local_sims  = local_sims)
  }
  # Adjacency Matrix
  W<- as(weight, "symmetricMatrix")
  W<- as.matrix(W/rowSums(W))

  m <- moran_I(x, y, W)
  # Global Moral
  global_moran <- m[[1]][1]
  # Local values
  m_i <- m[[2]]
  # local simulations
  sim<-simula_moran(x, y, W)
  local_sims <- sim$local_sims
  # global pseudo p-value
  # get all simulated global moran
  global_sims <- sim$global_sims

  # Proportion of simulated global values that are higher (in absolute terms) than the actual index
  moran_pvalue <- sum(abs(global_sims) > abs( global_moran )) / length(global_sims)

  # Identifying the significant values
  alpha2 <- .05  # for a 95% confidence interval
  probs <- c(alpha2/2, 1-alpha2/2)
  intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
  sig       <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )


  # Identifying the LISA clusters
  xp <- scale(x)[,1]
  yp <- scale(y)[,1]

  patterns <- as.character( interaction(xp > 0, W%*%yp > 0) )
  patterns <- patterns %>% str_replace_all("TRUE","High") %>% str_replace_all("FALSE","Low")

  patterns[sig==0] <- "Not significant"
  # map_sf$patterns <- patterns


  # Rename LISA clusters
  patterns2 <- factor(patterns, levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                      labels=c("High antigen - High tcell", "High antigen - Low tcell", "Low antigen - High tcell","Low antigen - Low tcell", "Not significant"))





  ### PLOT
  df<-data.frame(a=coords[,1],b=coords[,2])
  s<-ggplot(df, aes(b, -a)) +
      geom_point(aes(color=patterns2), size=3) +
      scale_color_manual(values=c("palevioletred2", "lightpink", "lightblue", "steelblue2", "grey80"))+theme_bw() +
      theme(axis.line = element_line(colour = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  print(s)
  png(paste0("plot/", sample_name, "_Spatial_Segmentation.png"))
  print(s)
  dev.off()






}
