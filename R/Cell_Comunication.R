
#' Cell_Communication
#'
#' @param pathlist a list containing the paths to CellphoneDB outdir
#' @param rm_complex logical value, infer whether delete the complex interaction in the final communication
#' @param top number, the top n ligand-receptor interaction in each cluster to be presented.
#'
#' @return a matrix, containing the ligand-receptor enrichment value
#' @export
#'
#'
Cell_Communication<-function(pathlist,rm_complex=FALSE){
  ###path refer to directory of pvalues and means
  all.interactionfile<-lapply(pathlist,Generate_Interactionfile)
  if(length(all.interactionfile) == 1){
    mat<-all.interactionfile[[1]]
  }else{
    for(i in 1:length(all.interactionfile)){
      if(i == 1){
        mat<-all.interactionfile[[i]]
      }else{
        mat<-rbind(mat,all.interactionfile[[i]])
      }
    }
  }
  #####generate data for plot
  if(rm_complex == TRUE){
    mat<-mat[-grep(pattern="complex",mat$interacting_pair),]
  }
  return(mat)
}


