#### this function is to merge pvalue and mean file in a data frame

Generate_Interactionfile<-function(path){
  pvalue<-read.delim(paste(path,"/pvalues.txt",sep = ""))
  mean<-read.delim(paste(path,"/means.txt",sep = ""))
  pvalue<-pvalue[,c(1,2,13,14)]
  group1<-colnames(pvalue)[3:4]
  pvalue<-tidyr::gather(pvalue,key = "cell-pair",value = "pvalue",group1)
  mean<-mean[,c(1,2,13,14)]
  group2<-colnames(mean)[3:4]
  mean<-tidyr::gather(mean,key = "cell-pair",value = "mean",group2)
  if(identical(group1,group2) == FALSE){
    stop("pvalue and mean file are not match.")
  }else{
    rownames(pvalue)<-paste(pvalue$id_cp_interaction,pvalue$`cell-pair`,sep="_")
    rownames(mean)<-paste(mean$id_cp_interaction,mean$`cell-pair`,sep="_")
    ###add mean information to pvalue dataframe
    pvalue$mean<-mean[rownames(pvalue),"mean"]
    pvalue<-pvalue[which(pvalue$pvalue <= 0.05),]
    pvalue<-pvalue[order(pvalue$pvalue),]
  }
  return(pvalue)
}

#### extract Target genes in a given sender and receiver celltypes
#' Extract_Target
#'
#' @param x a SpaTalk object
#' @param celltype_sender name of sender cell
#' @param celltype_receiver name of receiver cell
#' @param top top n ligand-receptor pairs to be extract
#' @import Seurat ggplot2
#'
#'
#' @return a series of characters containing targeted genes in the given sender and receiver cells
#' @export
#'
Extract_Target<-function(x,celltype_sender,celltype_receiver,top=20){
  ### extract LR pairs in the given sender and receiver cells
  lrpair<-x@lrpair
  lrpair<-lrpair[which(lrpair$celltype_sender == celltype_sender),]
  lrpair<-lrpair[which(lrpair$celltype_receiver == celltype_receiver),]
  lrpair <- lrpair[order(-lrpair$score), ]
  lrpair_target<-lrpair[1:top,1:2]
  ###extract LR target enriched path
  list_path<-list()
  for(i in 1:length(lrpair_target$ligand)){
    list_path[[paste(lrpair_target$ligand[i],lrpair_target$receptor[i],sep = "_")]]<-get_lr_path(object = x,celltype_sender = celltype_sender,celltype_receiver = celltype_receiver,ligand = lrpair_target$ligand[i],receptor = lrpair_target$receptor[i],min_gene_num = 0)
  }
  ### extract target genes
  list_pathgene<-c()
  for(i in names(list_path)){
    if(length(list_pathgene) == 0){
      list_pathgene<-unlist(list_path[[i]]$path_pvalue$gene)
    }else{
      list_pathgene<-c(list_pathgene,unlist(list_path[[i]]$path_pvalue$gene))
    }
  }
  ### modify gene format
  gene<-c()
  for(i in 1:length(names(list_path))){
    if(length(gene) == 0){
      gene<-stringr::str_split_fixed(list_pathgene,",",n=Inf)[i,]
    }else{
      gene<-c(gene,stringr::str_split_fixed(list_pathgene,",",n=Inf)[i,])
    }
  }
  return(gene)
}


### this function is used for GO/KEGG analysis and plot
Fun_Analysis<-function(gene,species,showCategory){
  if(species == "human"){
    OrgDb = "org.Hs.eg.db"
    organism="hsa"
  }else if (species == "mouse"){
    OrgDb = "org.Mm.eg.db"
    organism="mmu"
  }
  if(is.data.frame(gene)){
    GO <- clusterProfiler::compareCluster(gene~cluster,
                                          data = gene,
                                          keyType = "SYMBOL",
                                          fun = "enrichGO",
                                          OrgDb         = OrgDb,
                                          ont           = "BP",
                                          pAdjustMethod = "fdr",
                                          pvalueCutoff  = 0.2)
    entrez<-bitr(gene$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = OrgDb,drop = TRUE)
    entrez<-entrez[!duplicated(entrez$SYMBOL),]
    if(!is.null(GO)){
      p1<-dotplot(GO,showCategory = showCategory)+scale_fill_manual(values=paletteer_c("ggthemes::Classic Orange", 30))+ggtitle("GO")
      p1<-NULL
    }else{
      print("No Go term has been enriched in this gene set!")
    }
    gene$entrez<-NA
    for (num in 1:nrow(entrez)){
      gene[which(gene$gene == entrez$SYMBOL[num]),"entrez"] <-entrez$ENTREZID[num]
    }
    KEGG <- clusterProfiler::compareCluster(
      ENTREZID~cluster,
      data=gene,
      fun="enrichKEGG",
      organism=organism,
      pAdjustMethod = "fdr",
      pvalueCutoff=0.05)
    if(!is.null(KEGG)){
      p2<-dotplot(KEGG,showCategory = showCategory)+scale_fill_manual(values=paletteer_c("ggthemes::Classic Orange", 30))+ggtitle("KEGG")
    }else{
      print("No KEGG term has been enriched in this gene set!")
      p2<-NULL
    }
  }else if(is.character(gene)){
    GO <- clusterProfiler::enrichGO(gene,
                                    OrgDb = OrgDb,
                                    keyType = "SYMBOL",
                                    ont="BP"
                                    )
    if(!is.null(GO)){
      p1<-barplot(GO,showCategory=showCategory)+scale_fill_gsea(reverse = T)+ggtitle("GO")
    }else{
      print("No GO term has been enriched in this gene set!")
      p1<-NULL
    }
    entrez<-bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = OrgDb,drop = TRUE)
    KEGG <-clusterProfiler::enrichKEGG(unique(entrez$ENTREZID),
                                       organism = organism,
                                       keyType = "kegg"
    )
    if(!is.null(KEGG)){
      p2<-barplot(KEGG,showCategory=showCategory)+scale_fill_gsea(reverse = T)+ggtitle("KEGG")
    }else{
      print("No KEGG term has been enriched in this gene set!")
      p2<-NULL
    }
  }
  p<-cowplot::plot_grid(plotlist = list(p1,p2))
  print(p)
  return(list(G0=GO,KEGG=KEGG))
}

####
Qc_Check<-function(x,outdir,spatial=FALSE){
  x <- NormalizeData(x)
  x <-FindVariableFeatures(x,nfeatures=2000)
  x <- ScaleData(object = x,assay = "RNA",features=VariableFeatures(x))
  x <- RunPCA(object = x,assay = "RNA",features=VariableFeatures(x))
  x <- FindNeighbors(x, reduction = "pca")
  ####vln plot
  pdf(paste(outdir, "Vlnplot.pdf", sep = "/"), width = 6, height = 4)
  p <- VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, combine = F)
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoLegend() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45))
  }
  p <- cowplot::plot_grid(plotlist = p, ncol = 3)
  print(p)
  dev.off()
  #####feature plot
  pdf(paste(outdir, "umap_featureplot.pdf", sep = "/"), width = 14, height = 6)
  p <- FeaturePlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), combine = F)
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45))
  }
  print(cowplot::plot_grid(plotlist = p, ncol = 3))
  dev.off()
  ##### spatial feature plot
  if(spatial == FALSE){
    if("spatial" %in% names(x@reductions)){
      pdf(paste(outdir, "spatial_featureplot.pdf", sep = "/"), width = 14, height = 6)
      p <- FeaturePlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), combine = F,reduction = "spatial")
      for (i in 1:length(p)) {
        p[[i]] <- p[[i]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45))
      }
      print(cowplot::plot_grid(plotlist = p, ncol = 3))
      dev.off()
    }
  }else{
    pdf(paste(outdir, "spatial_featurplot.pdf", sep = "/"), width = 7, height = 7)
    p <- SpatialFeaturePlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), combine = F)
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45))
    }
    print(cowplot::plot_grid(plotlist = p, ncol = 3))
    dev.off()
  }
}
###
