#' Standard_Spatalk
#'
#' @param sp a seurat object or count matrix  of spatial transcriptome data
#' @param sc a seurat object or count matrix of reference scRNAseq data
#' @param sc_celltype a series character label of the sc reference's cell barcodes
#' @param sp_meta a matrix recording the spatial locate information, column as x, y and rownames are spot barcodes
#' @param geneinfo the same as description in SpaTalk
#' @param lrpairs the same as description in SpaTalk
#' @param pathways the same as description in SpaTalk
#'
#' @return a spatalk object with celltype and ligand-receptor enrichment information
#' @export
#'

Standard_Spatalk<-function(sp,sc,sc_celltype,sp_meta=NULL,geneinfo,lrpairs,pathways){
  if(!is(sp,"Seurat") & !is(sp,"matrix")){
    stop("sp should be a seurat object or count matrix")
  }else if (!is(sc,"Seurat") & !is(sc,"matrix")) {
    stop("sc should be a seurat object or count matrix")
  }
  ##### retrieve sp count and spatial location
  if(is(sp,"Seurat")){
   st_data<-sp@assays$RNA@counts
   #### if sp is from stereoseq
   st_meta<-sp@reductions$spatial@cell.embeddings
   st_meta<-as.data.frame(st_meta)
   st_meta$spot<-rownames(st_meta)
   st_meta<-st_meta[,c(3,1,2)]
   colnames(st_meta)[c(2,3)]<-c("x","y")
  }
  else{
    if(is.null(sp_meta)){
      stop("sp_meta is required!")
    }else{
      st_data<-sp
      st_meta<-sp_meta
    }
  }
  #####retrieve scRNAseq reference count and cell annotation
  if(is(sc,"Seurat")){
    sc_data<-sc@assays$RNA@counts
  }else{
    sc_data<-sc
  }
  st_data <- rev_gene(data = as.matrix(st_data),
                       data_type = "count",
                       species = "Human",
                       geneinfo = geneinfo)
  obj <- createSpaTalk(st_data = as.matrix(st_data),
                       st_meta = st_meta,
                       species = "Human",
                       if_st_is_sc = F,
                       spot_max_cell = 10)
  ######deconvolution using spatalk default methods
  obj <- dec_celltype(object = obj,
                      sc_data = as.matrix(sc_data),
                      sc_celltype = sc_celltype,
                      use_n_cores=20,
                      method = 1)
  ##### find ligand receptor pair
  obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
  obj <- dec_cci_all(object)
}

