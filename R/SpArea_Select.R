#' SpArea_Select
#' this function is to select boundary spot manually.
#' @param x a seurat object
#' @param ident ident to label the subset object (eg. "xx_boundary")
#' @param group group to cluster your data, should be in your meta.data (eg. celltype slot in your seurat meta.data)
#' @param interestsID characters infer your interested spot to select, the character should be in the group above.
#' @param cols Vector of colors, each color corresponds to an identity class.
#' @param spatial10x logical value to infer if the data is from cellranger output
#' @param saveCPDB_dir directory path to save the files for CellphoneDB
#'
#' @return a seurat object within selected spots, present the selected area in the plot
#' @export
#'
SpArea_Select<-function(x,ident,group,interestsID,cols=NULL,spatial10x=FALSE,saveCPDB_dir=NULL){
  if(group %in% colnames(x@meta.data) == FALSE){
    stop(paste(group,"does not exist in your meta.data",sep=" "))
  }else{
    if(all(interestsID %in% x@meta.data[[group]]) == FALSE){
      outid<-setdiff(interestsID,x@meta.data[[group]])
      stop(paste(outid,"is not in your selected group!",sep = " "))
    }
  }
  if(spatial10x == FALSE){
    p<-DimPlot(x,group.by=group,reduction = "spatial",cols=cols)
  }else{
    p<-SpatialDimPlot(x,group.by=group,cols=cols)
  }
  select_X<-CellSelector(p,object = x,ident = ident)
  select_X<-subset(select_X,idents= ident)
  Idents(select_X)<-select_X@meta.data[[group]]
  outlier<-which(levels(Idents(select_X)) %in%  interestsID== FALSE)
  levels(Idents(select_X))[outlier]<-NA
  ###recheck spot to remove spot not in boundary region
  mannual=T
  while(mannual == T){
    if(spatial10x == FALSE){
      p<-DimPlot(select_X,reduction = "spatial")
    }else{
      p<-SpatialDimPlot(select_X)
    }
    select_X<-CellSelector(p,object = select_X,ident = NA)
    mannual<-readline(prompt="Keep select outlier spot? (T or F).  ")
  }
  ### plot the final selected area
  if(spatial10x == FALSE){
    p<-DimPlot(select_X,reduction = "spatial")+theme_void()
  }else{
    p<-SpatialDimPlot(select_X)+theme_void()
  }
  print(p)
  if(!is.null(saveCPDB_dir)){
    select_X<-subset(select_X,idents=levels(Idents(select)))
    select_X$cell<-rownames(select_X@meta.data)
    select_X$celltype<-Idents(select_X)
    df = select_X@meta.data[, c('cell', 'celltype')]
    ## merge interested id together as a new interact id
    name=paste0(levels(Idents(select_X)),collapse = "_")
    ## in stereoseq data the barcode ID are numbers  a "X" is needed in front each number
    df$cell<-paste("X",df$cell,sep = "")
    write.table(df,paste(saveCPDB_dir,"/",name,"_meta.txt",sep = ""),sep = '\t', quote = F, row.names = F)
    write.table(select_X@assays$RNA@data,paste(saveCPDB_dir,"/",name,"_count.txt",sep = ""),sep = "\t",quote = F)
  }
  return(select_X)
}
