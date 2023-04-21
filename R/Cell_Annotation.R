#' Title
#'
#' @param x a seurat object
#' @param method a character, only can set to "singleR" in this version
#' @param ref reference single cell RNAseq data
#' @param label reference labels
#' @param cluster logical value, whether to label cells by cluster
#'
#' @return a seurat object with predicted cell label
#' @import singleR, Seurat
#' @export
#'
Cell_Annotation<-function(x,method="singleR",ref,label,cluster=TRUE){
  if(method == "singleR"){
    test_for_singler<-GetAssayData(x, slot="data")
    if(cluster == FALSE){
      singler1<-SingleR(test= test_for_singler , ref = ref ,labels = label)
      x$singleR_predicion<-x$labels
    }
    else if(cluster == TRUE){
      res=colnames(x@meta.data)[grep(colnames(x@meta.data),pattern = "RNA_snn_res")]
      for(i in res){
      clusters=x@meta.data[[i]]
      names(clusters)<-rownames(x@meta.data)
      singler2<-SingleR(test= test_for_singler,ref = ref,labels = label,clusters=clusters)
      celltype = data.frame(ClusterID=rownames(singler2), celltype=singler2$labels, stringsAsFactors = FALSE)
      for(num in 1:nrow(celltype)){
        ##### names clusters by singleR predicted ID
        x@meta.data[which(clusters == celltype$ClusterID[num]),paste(i,"singleR_predictions",sep="_")] <- celltype$celltype[num]
        }
      }
    }
  }else if(method == "Seurat"){
    anchors <- Seurat::FindTransferAnchors(reference = ref, query = x, verbose = F)
    celltype.predictions<- Seurat::TransferData(anchorset = anchors, refdata = label,
                                                weight.reduction = x[["pca"]],dims = 1:30)
  }else{
    print("this method is under development in our script")
  }
  return(x)
}
