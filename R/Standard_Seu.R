#' Standard_Seu
#'
#' @param x a seurat object, a count matrix or a parent directory path to cellranger filtered_feature_bc_matrix
#' @param spatial_dir the cellranger out spatial directory. Due to lack of spatial HE images, if your data is from stereoseq, this parameter should set to NULL.
#' @param QC_dir if not null, should be the QC outdir path you want save your QC figures. Make sure you have access to the QC directory
#' @param project your sample name for identification
#' @param filter logical, whether to filter cells/spots in top5% in lowest minfeature and highest maxfeature with top5% mitochondrial genes.
#'
#' @return a seurat object
#' @import Seurat ggplot2
#' @export
#'

Standard_Seu<-function(x,spatial_dir=NULL,QC_dir=NULL,project,filter=TRUE){
  if(is.matrix(x) || is(x, "dgCMatrix")){
    x<-Seurat::CreateSeuratObject(counts = x,min.cells=0,min.features=0,project = project)
  }else if(is.character(x)){
    x<-Seurat::Read10X(data.dir = paste(x, "filtered_feature_bc_matrix", sep = "/"))
    x<-Seurat::CreateSeuratObject(counts = x,min.cells=0,min.features=0,project = project)
  }else if(is(x,"Seurat")){
    x<-x
  }else{
    print(paste(typeof(x), "is not supported. your input should be a count matrix, a 10x matrix directory path or a seurat object!"))
  }
  if(!is.null(spatial_dir)){
    Ximage <- Seurat::Read10X_Image(image.dir = spatial_dir)
    Seurat::DefaultAssay(Ximage) <- "RNA"
    # link matrix and image file
    Ximage <- Ximage[colnames(x)]
    x[["image"]] <- Ximage
    spatial=TRUE
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  }else{
    spatial=FALSE
    #### if "total_counts" in x meta.data, means that x is a data from stereotype
    if("total_counts" %in% colnames(x@meta.data)){
      colnames(x@meta.data)[which(colnames(x@meta.data) == "total_counts")]<-"nCount_RNA"
      colnames(x@meta.data)[which(colnames(x@meta.data) == "pct_counts_mt")]<-"percent.mt"
      colnames(x@meta.data)[which(colnames(x@meta.data) == "n_genes_by_counts")]<-"nFeature_RNA"
    }
  }
  if(!is.null(QC_dir)){
    if(!is.character(QC_dir)){
      stop("QC_dir should be the outdir path for Quality control results")
    }else{
      Qoutdir=paste(QC_dir,"QC_dir",sep = "/")
      if(!dir.exists(Qoutdir)){
        dir.create(Qoutdir)
      }
      outdir=paste(Qoutdir,project,sep = "/")
      if(!dir.exists(outdir)){
        dir.create(outdir)
      }
      Qc_Check(x,outdir,spatial=spatial)
    }
  }
  if(filter == TRUE){
    minFeature<-quantile(x$nFeature_RNA,probs = 0.05)
    maxFeature<-quantile(x$nFeature_RNA,probs = 0.95)
    maxMt<-quantile(x$percent.mt,probs = 0.95)
    x <- subset(x, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & percent.mt < maxMt)
  }
  x <- NormalizeData(x)
  x <-FindVariableFeatures(x,nfeatures=2000)
  x <- ScaleData(object = x,assay = "RNA",features=VariableFeatures(x))
  x <- RunPCA(object = x,assay = "RNA",features=VariableFeatures(x))
  x <- FindNeighbors(x, reduction = "pca")
  for(i in seq(0.1,0.9,by=0.1)){
    x<-FindClusters(x,algorithm = 1,resolution = i)
  }
  x<-RunUMAP(x,dims=1:20)
  x<-RunTSNE(x,check_duplicates = FALSE)
  return(x)
}
