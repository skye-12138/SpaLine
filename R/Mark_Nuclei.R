
#' Mark_Nuclei
#'
#' @param path a character infer the parent directory to cellranger filtered_feature_bc_matrix
#' @param QC_dir a character indicate the path you want save your QC figures. Make sure you have access to the QC directory
#' @param project a character indicate your sample name for identification
#'
#' @return a seurat object within nuclei information in meta.data
#' @export
#'

Mark_Nuclei<-function(path,QC_dir=NULL,project){
#  env<-reticulate::conda_list()
#  if("SpaLine_conda" %in% env$name){
#    reticulate::use_condaenv("SpaLine_conda", required = TRUE)
#  }else{
#    reticulate::conda_create("SpaLine_conda")
#    reticulate::use_condaenv("SpaLine_conda", required = TRUE)
#  }
#  reticulate::py_config()
  reticulate::py_run_file(system.file("python","st_nucleus.py",package = "SpaLine"))
  reticulate::py_run_file(system.file("python","st_segment.py",package = "SpaLine"))
  adata <- feature_tiling(inDir=path,outDir=QC_dir,project=project)
  ###feature extract
  feature_mtx <- adata$obsm['X_tile_feature']
  obs_names <- rownames(adata$obs)
  rownames(feature_mtx) <- obs_names
  colnames(feature_mtx) <- paste("feature",seq_len(ncol(feature_mtx)),sep = "")
  ###Creat Seu.obj
  x<-Standard_Seu(t(feature_mtx), project = project,spatial_dir=paste(path,"spatial",sep = "/"))
  #######extract nuclei information
  adataF <- morph_watershed(adata,copy = T)
  n_nuclei <- adataF$obs[,"n_nuclei"]
  names(n_nuclei) <- rownames(adataF$obs)
  x@meta.data[,"n_nuclei"] <- n_nuclei[match(rownames(x@meta.data),names(n_nuclei))]
  #########
  if(!is.null(QC_dir)){
    #### different resolution for umap plot
    png(paste(QC_dir,"Morph Different resolution.png",sep = "/"),width=3500,height=2100,res=200)
    DiffResolution_plot <- lapply(seq(0.1,1,by=0.1),function(i){
      p <- SpatialDimPlot(x,label=F,group.by = paste("RNA_snn_res.",i,sep=""),cols = .cluster_cols)+
        labs(title=paste("Resolution = ",i,sep=""))
      return(p)
    })
    p_res <- cowplot::plot_grid(plotlist = DiffResolution_plot,ncol=5)
    print(p_res)
    dev.off()
    #### nuclei number in each cluster under different resolution
    png(paste(QC_dir,"nuclei violin Different resolution.png",sep = "/"),width=6000,height=6000,res=200)
    DiffResolution_plot <- lapply(seq(0.1,1,by=0.1),function(i){
      p <- VlnPlot(x,features = c("n_nuclei"),group.by = paste("RNA_snn_res.",i,sep=""),pt.size = 1,cols = .cluster_cols)+NoLegend()+
        labs(title=paste("Resolution = ",i,sep=""))
      return(p)
    })
    p_res <- cowplot::plot_grid(plotlist = DiffResolution_plot,ncol=1)
    print(p_res)
    dev.off()
    ####nuclei distribution in spatial spot
    png(paste(QC_dir,"nuclei in spatial.png",sep = "/"))
    p<-SpatialFeaturePlot(x,features = "n_nuclei",pt.size.factor = 1.2,alpha = c(0,1))
    print(p)
    dev.off()
  }
  return(x)
}
