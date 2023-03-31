
Feature_plot<-function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
  c("lightgrey", "#ff0000", "#00ff00")
} else {
  c("#EBDCDD","#AE123A")
}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA,
reduction = NULL, split.by = NULL, keep.scale = "feature",
shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5,
label = FALSE, label.size = 4, label.color = "black", repel = FALSE,
ncol = NULL, coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL,
interactive = FALSE, combine = TRUE, raster = NULL, raster.dpi = c(512,512)) {
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ",
            "parameter instead for equivalent functionality.",
            call. = FALSE, immediate. = TRUE)
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(object = object, feature = features[1],
                        dims = dims, reduction = reduction, slot = slot))
  }
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature",
                                                        "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
                                                                                           size = 14, margin = margin(r = 7)))
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
                   `0` = {
                     warning("No colors provided, using default colors",
                             call. = FALSE, immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
                             call. = FALSE, immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '",
                             default.colors[1], "' for double-negatives",
                             call. = FALSE, immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three",
                             call. = FALSE, immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(dims, "ident",
                                              features), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", paste(features,
                                                              collapse = ", "), " in slot ", slot, call. = FALSE)
  }
  else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[,
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[,
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
  ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data),
                                     FUN = function(index) {
                                       data.feature <- as.vector(x = data[, index])
                                       min.use <- SetQuantile(cutoff = min.cutoff[index -
                                                                                    3], data.feature)
                                       max.use <- SetQuantile(cutoff = max.cutoff[index -
                                                                                    3], data.feature)
                                       data.feature[data.feature < min.use] <- min.use
                                       data.feature[data.feature > max.use] <- max.use
                                       if (brewer.gran == 2) {
                                         return(data.feature)
                                       }
                                       data.cut <- if (all(data.feature == 0)) {
                                         0
                                       }
                                       else {
                                         as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                                                                          breaks = brewer.gran)))
                                       }
                                       return(data.cut)
                                     })
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  }
  else {
    switch(EXPR = split.by, ident = Idents(object = object)[cells,
                                                            drop = TRUE], object[[split.by, drop = TRUE]][cells,
                                                                                                          drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[,
                                                                   dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[,
                                                               dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold,
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ],
                   as.vector(x = color.matrix))
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident,
                      , drop = FALSE]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[,
                                                       features]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ",
             paste(no.expression, collapse = ", "), call. = FALSE)
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")],
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[,
                                                              feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      }
      else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", feature,
                                   shape.by)]
      plot <- SingleDimPlot(data = data.single, dims = dims,
                            col.by = feature, order = order, pt.size = pt.size,
                            cols = cols.use, shape.by = shape.by, label = FALSE,
                            raster = raster, raster.dpi = raster.dpi) + scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) + theme_cowplot() +
        CenterTitle()
      if (label) {
        plot <- LabelClusters(plot = plot, id = "ident",
                              repel = repel, size = label.size, color = label.color)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA,
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        }
        else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident),
                                                                    limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(),
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(),
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                               axis.title.x = element_blank())
        }
      }
      else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        }
        else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (",
                    unique.feature.exp, ") of ", feature, ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            }
            else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad,
                                                                       guide = "colorbar"))
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" &&
          !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols,
                                                              limits = c(min.feature.value, max.feature.value)))
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots,
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                                                                                                                                  1, yes = levels(x = data$split)[ii], no = "")),
                                                                                              expand = c(0, 0)) + labs(x = features[1], y = features[2],
                                                                                                                       title = if (ii == 1) {
                                                                                                                         paste("Color threshold:", blend.threshold)
                                                                                                                       } else {
                                                                                                                         NULL
                                                                                                                       }) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol,
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  }
  else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""),
                                                                   limits = ylims) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]),
                                                            limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) ==
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      }
      else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots),
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      #### add theme in each plot --differ from seurat
      plots <-lapply(plots, function(x){x+theme_void()})
      plots <- patchwork::wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    }
    else {
      #### add theme in each plot --differ from seurat
      plots <-lapply(plots, function(x){x+theme_void()})
      plots <- patchwork::wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
                                       length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" &&
        !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(colors = cols,
                                                              limits = c(min.feature.value, max.feature.value)))
    }
  }
  return(plots)
}



#' CellRatio_bar
#'
#' @param x spatalk object from Standard_Spatalk or a seurat object
#' (this function will calculate the idents of the seurat object's ratio)
#'
#' @export
#'

CellRatio_bar<-function(x,meta=NULL){
  if(is(x,"SpaTalk")){
    mat<-as.data.frame(table(x@meta$newmeta$celltype))
  }else if(is(x,"Seurat")){
    if(is.null(meta)){
      mat<-as.data.frame(table(Idents(x)))
    }else{
      mat<-as.data.frame(table(x@meta.data[[meta]]))
    }
  }

  ### calculate cell ratio
  mat$ratio<-mat$Freq/sum(mat$Freq)
  mat$Var1<-factor(mat$Var1,levels =mat$Var1[order(mat$ratio)])
  p<-ggplot(mat) +
    aes(x = Var1, fill = Var1, weight = ratio) +
    geom_bar() +
    scale_fill_hue(direction = 1) +
    theme_minimal()
  #+scale_fill_manual(
      #values = color.cell[levels(celltype_ratio_p160$Var1)])
  print(p)
  return(p)
}




#' LR_plots
#'
#' @param mat matrix get from Cell_Comunication
#' @param mat matrix get from Cell_Comunication
#'
#' @export

LR_plots<-function(mat,weight="pvalue"){
  if(weight == "pvalue"){
    mat<-mat %>% group_by(`cell-pair`) %>% top_n(n=-top,wt=pvalue)
  }else if(weight == "mean"){
    mat<-mat %>% group_by(`cell-pair`) %>% top_n(n=-top,wt=mean)
  }
  mat$pvalue[which(mat$pvalue <= 0.05 & mat$pvalue >0.01)]<-1
  mat$pvalue[which(mat$pvalue <= 0.01 & mat$pvalue >0.001)]<-2
  mat$pvalue[which(mat$pvalue <= 0.001)]<-3
  mat$pvalue<-as.factor(mat$pvalue)
  return(mat)
  p<-ggplot(mat,aes(x=mat$`cell-pair`,y=mat$interacting_pair))+
    geom_point(aes(size=mat$pvalue,
                   color=mat$mean))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
    scale_size_discrete(range = c(1,2,3),labels=c("*","**","***"))+
    labs(x=NULL,y=NULL)+scale_colour_gradientn(colours = paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30,direction = -1))+labs(color="mean",size="significant")
  print(p)
}


#' LRpair_enrich
#'
#' @param object SpaTalk object
#' @param celltype_sender name of sender cell
#' @param celltype_receiver name of receiver cell
#' @param top_lrpairs number of top lrpairs for plotting
#' @param color color for the cells in heatmap.
#' @param border_color color of cell borders on heatmap, use NA if no border should be drawn.
#' @param type recommend set to "number", plot significant LR with number as fill color
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param number_color color of the text
#' @param color_low 	for 'number' type, define the color for the lowest value.
#' @param color_high for 'number' type, define the color for the highest value
#'
#' @export
#'
LRpair_enrich<-function (object, celltype_sender, celltype_receiver, top_lrpairs = 20,
                         color = NULL, border_color = "black", type = "number", fontsize_number = 5,
                         number_color = "black", color_low = NULL, color_high = NULL)
{
  if (!is(object, "SpaTalk")) {
    stop("Invalid class for object: must be 'SpaTalk'!")
  }
  st_type <- object@para$st_type
  if_skip_dec_celltype <- object@para$if_skip_dec_celltype
  if (st_type == "single-cell") {
    st_meta <- object@meta$rawmeta
    if (if_skip_dec_celltype) {
      st_data <- object@data$rawdata
    }
    else {
      st_data <- object@data$rawndata
    }
  }
  else {
    if (if_skip_dec_celltype) {
      st_meta <- object@meta$rawmeta
      colnames(st_meta)[1] <- "cell"
      st_data <- object@data$rawdata
    }
    else {
      st_meta <- object@meta$newmeta
      st_data <- object@data$newdata
    }
  }
  if (!celltype_sender %in% st_meta$celltype) {
    stop("Please provide the correct name of celltype_sender!")
  }
  if (!celltype_receiver %in% st_meta$celltype) {
    stop("Please provide the correct name of celltype_receiver!")
  }
  if (is.null(color)) {
    cellname <- unique(st_meta$celltype)
    cellname <- cellname[order(cellname)]
    if ("unsure" %in% cellname) {
      cellname <- cellname[-which(cellname == "unsure")]
    }
    col_manual <- ggpubr::get_palette(palette = "lancet",
                                      k = length(cellname))
    color <- col_manual[which(cellname == celltype_receiver)]
  }
  heat_col <- (grDevices::colorRampPalette(c("white", color)))(2)
  cell_pair <- object@cellpair
  cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
  if (is.null(cell_pair)) {
    stop("No LR pairs found from the celltype_sender to celltype_receiver!")
  }
  cell_pair <- cell_pair[cell_pair$cell_sender %in% st_meta$cell &
                           cell_pair$cell_receiver %in% st_meta$cell, ]
  lrpair <- object@lrpair
  lrpair <- lrpair[lrpair$celltype_sender == celltype_sender &
                     lrpair$celltype_receiver == celltype_receiver, ]
  lrpair <- lrpair[order(-lrpair$score), ]
  if (nrow(lrpair) > top_lrpairs) {
    lrpair <- lrpair[1:top_lrpairs, ]
  }
  ligand <- unique(lrpair$ligand)
  receptor <- unique(lrpair$receptor)
  if (type == "sig") {
    lrpair_real <- object@lr_path$lrpairs
    lrpair_real <- lrpair_real[, c(1, 2)]
    lrpair_real$score <- 1
    lrpair_mat <- reshape2::dcast(lrpair_real, formula = ligand ~
                                    receptor, fill = 0, value.var = "score")
    rownames(lrpair_mat) <- lrpair_mat$ligand
    lrpair_mat <- lrpair_mat[, -1]
    lrpair_mat <- lrpair_mat[ligand, receptor]
    if (!is.data.frame(lrpair_mat)) {
      stop("Limited number of ligand-receptor interactions!")
    }
    plot_res <- matrix("", nrow = length(ligand), ncol = length(receptor))
    rownames(plot_res) <- ligand
    colnames(plot_res) <- receptor
    for (i in 1:nrow(lrpair)) {
      plot_res[lrpair$ligand[i], lrpair$receptor[i]] <- "*"
    }
    pheatmap::pheatmap(lrpair_mat, cluster_cols = F, cluster_rows = F,
                       color = heat_col, border_color = border_color, legend = F,
                       display_numbers = plot_res, fontsize_number = fontsize_number,
                       number_color = number_color, main = "Significantly enriched LRI pairs")
  }
  else {
    if (is.null(color_low)) {
      color_low <- "orange"
    }
    if (is.null(color_high)) {
      color_high <- "red"
    }
    lrpair_real <- lrpair[, c("ligand", "receptor", "lr_co_exp_num")]
    lrpair_mat <- reshape2::dcast(lrpair_real, formula = ligand ~
                                    receptor, fill = 0, value.var = "lr_co_exp_num")
    rownames(lrpair_mat) <- lrpair_mat$ligand
    lrpair_mat <- lrpair_mat[, -1]
    lrpair_mat <- lrpair_mat[ligand, receptor]
    plot_res <- matrix("", nrow = length(ligand), ncol = length(receptor))
    rownames(plot_res) <- ligand
    colnames(plot_res) <- receptor
    for (i in 1:nrow(lrpair)) {
      plot_res[lrpair$ligand[i], lrpair$receptor[i]] <- "*"
    }
    heat_color <- (grDevices::colorRampPalette(c(color_low,
                                                 color_high)))(max(as.matrix(lrpair_mat)) - 1)
    heat_color <- c("white", heat_color)
    pheatmap::pheatmap(lrpair_mat, cluster_cols = F, cluster_rows = F,
                       border_color = border_color, color = heat_color,
                       display_numbers = plot_res,
                       fontsize_number = fontsize_number,
                       number_color = number_color,
                       annotation_names_row=TRUE,
                       annotation_names_col=TRUE,
                       main = "Number of spatial LRI pairs")
  }
}

#' Fun_Plot
#'
#' @param x a SpaTalk object or data.frame with gene and cluster information or a series of genes
#' @param species name of species
#' @param showCategory number of GO/KEGG categories to be shown
#' @param celltype_sender name of sender cell
#' @param celltype_receiver name of receiver cell
#' @param top top n ligand-receptor pairs to be extract
#'
#' @return a list containing GO and KEGG output
#' @export
#'

Fun_Plot<-function(x,species,showCategory=20,celltype_sender=NULL,celltype_receiver=NULL,top=20){
  if(is.character(x) || is.data.frame(x)){
    x<-x
  }else if(is(x,"SpaTalk")){
    if(is.null(celltype_sender) || is.null(celltype_receiver)){
      stop("sender and reveiver celltypes is needed !")
    }else{
      x<-Extract_Target(x,celltype_sender=celltype_sender,celltype_receiver=celltype_receiver,top=top)
    }
  }else{
    stop("the input is not in a proper format!")
  }
  out<-Fun_Analysis(x,species=species,showCategory=showCategory)
  return(out)
}
