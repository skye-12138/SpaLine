#' Title
#'
#' @param x a seurat object
#' @param method a character, only can set to "SingleR" in this version
#' @param ref reference single cell RNAseq data
#' @param label reference labels
#' @param cluster logical value, whether to label cells by cluster
#'
#' @return a seurat object with predicted cell label
#' @import SingleR, Seurat
#' @export
#'
Cell_Annotation<-function(x,method="SingleR",ref,label,cluster=TRUE,marker=NULL,weight_id = "avg_log2FC",scale=TRUE){
  if(method == "SingleR"){
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
        x@meta.data[which(clusters == celltype$ClusterID[num]),paste(i,"SingleR_predictions",sep="_")] <- celltype$celltype[num]
        }
      }
    }
  }else if(method == "Seurat"){
    x<-.seurat(sc_data=ref,st_data=x,sc_celltype=label)
  }else if(method == "SPOTlight"){
    splist<-.spotlight(sc_data=ref,st_data=x,sc_celltype=label,marker=marker,weight_id = "avg_log2FC",scale=TRUE)
  }
  else{
    print("this method is under development in our script")
  }
  return(x)
}
