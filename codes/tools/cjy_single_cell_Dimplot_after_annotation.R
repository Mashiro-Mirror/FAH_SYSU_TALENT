cjy_single_cell_Dimplot_after_annotation = function(scRNA, annotation_column, 
                                                           ggtheme = theme_void() + 
                                                             theme(plot.margin = margin(5.5,15,5.5,5.5)),
                                                           pt.size = 0.4, linetype = "solid", font.size = 4, 
                                                           legend.point.size = 5,
                                                           legend.font.size = 12,legend.title.size = 14,
                                                    elipse_or_not = T, anno.label = T,
                                                    mycol = c4a("bold", length(unique(celltype)))){
  ## Require the necessary packages
  require(ggplot2)
  require(Seurat)
  require(tidyverse)
  require(cols4all)
  require(tidydr)                                                                
  ## Extract the umap coordinates from the scRNA object
  umap = as.data.frame(scRNA@reductions$umap@cell.embeddings)
  ## Extract the annotation column from the scRNA object
  celltype = scRNA@meta.data[, annotation_column]
  ## Combine the umap coordinates and the annotation column into a single data frame
  umap_celltype = cbind(umap, celltype)
  ## Prepare annotation label at the center of each cluster(using the mean of the umap coordinates)
  label = umap_celltype %>% 
    group_by(celltype) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  ## Prepare my own color:
  mycol = mycol
  #####plotting umap with cell type annotation using ggplot##########
  if (elipse_or_not) {
    p = ggplot(data=umap_celltype, aes(x=UMAP_1, y=UMAP_2)) + 
      ## First Add elipse confidence interval
      stat_ellipse(aes(color = celltype), level = 0.95, linetype = ifelse(linetype == "solid", 1, 2), show.legend = FALSE, geom = "polygon", alpha = .1) + 
      ## add the basic theme here
      ggtheme +
      ## Second add the points (Considering the relationship among the ggplot layers, the points should be added after the elipse)
      geom_point(aes(color = celltype),size=pt.size) + 
      ## Third add umap arrow
      theme_dr(xlength = 0.2, # axis x length
               ylength = 0.2, # axis y length
               arrow = grid::arrow(length = unit(.1, "inches"), ## arrow size/length
                                   ends = "last", type = "closed")) + ## arrow type
      theme(panel.grid = element_blank(), legend.text = element_text(size = legend.font.size),
            legend.title = element_text(size = legend.title.size)) + ## remove the grid
      ## Fourth set the size of the point in the legend
      guides(color = guide_legend(override.aes = list(size = legend.point.size),
                                  title = annotation_column)) +
      ## Fifth try to use my own color
      scale_color_manual(values = mycol) +
      scale_fill_manual(values = mycol)
  } else {
    p = ggplot(data=umap_celltype, aes(x=UMAP_1, y=UMAP_2)) + 
      ## First Add elipse confidence interval
      # stat_ellipse(aes(color = celltype), level = 0.95, linetype = ifelse(linetype == "solid", 1, 2), show.legend = FALSE, geom = "polygon", alpha = .1) + 
      ## add the basic theme here
      ggtheme +
      ## Second add the points (Considering the relationship among the ggplot layers, the points should be added after the elipse)
      geom_point(aes(color = celltype),size=pt.size) + 
      ## Third add umap arrow
      theme_dr(xlength = 0.2, # axis x length
               ylength = 0.2, # axis y length
               arrow = grid::arrow(length = unit(.1, "inches"), ## arrow size/length
                                   ends = "last", type = "closed")) + ## arrow type
      theme(panel.grid = element_blank(), legend.text = element_text(size = legend.font.size),
            legend.title = element_text(size = legend.title.size)) + ## remove the grid
      ## Fourth set the size of the point in the legend
      guides(color = guide_legend(override.aes = list(size = legend.point.size),
                                  title = annotation_column)) +
      ## Fifth try to use my own color
      scale_color_manual(values = mycol) +
      scale_fill_manual(values = mycol)
  }
  if (anno.label){
    p = p + geom_text(data = label, aes(x = UMAP_1, y = UMAP_2, label = celltype), size = font.size, color = "black", fontface = "bold")
  } else {
    p = p
  }
  return(p)
}
