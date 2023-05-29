
.seurat<-function(sc_data,st_data,sc_celltype){
  #sc_data<-SpaLine::Standard_Seu(x = sc_data,project = "sc",filter=FALSE)
  sc_data$celltype <- sc_celltype
  # Seurat intergration
  #st_data<- Standard_Seu(x = st_data,project = "st",filter=FALSE)
  anchors <- Seurat::FindTransferAnchors(reference = sc_data, query = st_data, verbose = F)
  celltype.predictions <- Seurat::TransferData(anchorset = anchors, refdata = sc_data$celltype,
                         weight.reduction = st_data[["pca"]],dims = 1:30,k.weight = 10)
  colnames(celltype.predictions)<-gsub(pattern = "predicted",replacement = "Seurat_predicted",x = colnames(celltype.predictions))
  colnames(celltype.predictions)<-gsub(pattern = "prediction",replacement = "Seurat_predicted",x = colnames(celltype.predictions))
  st_data<-Seurat::AddMetaData(st_data,metadata = celltype.predictions)
  return(st_data)
}

.cell2location<-function(){

}

.spotlight<-function(sc_data,st_data,sc_celltype,marker,weight_id = "avg_log2FC",scale=TRUE){
  x<-SPOTlight::SPOTlight(x=sc_data,y=st_data@assays$RNA@counts,groups = sc_celltype,mgs = marker,weight_id=weight_id,scale=scale)
  return(x)
}

.music<-function(){

}

.card<-function(){

}
