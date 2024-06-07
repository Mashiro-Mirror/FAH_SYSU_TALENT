cjy_single_cell_chisq_heatmap = function(single.cell.obj, column.to.chisq, compare.groups,
                                         heatmap_palette = colorRampPalette(c("#FFFFFF",'#FFCC66',"#CC3300"))(100), 
                                         heatmap_path = "chisq_heatmap/",
                                         cellwidth = 60, cellheight = 32, fontsize = 24, width = 20, height = 12, 
                                         display_pval_in_heatmap = FALSE, using_complexheatmap = TRUE) {
  require(Seurat)
  require(tidyverse)
  require(circlize)
  ## Make sure the heatmap_path ends with "/" and exists
  if (!endsWith(heatmap_path, "/")) {
    heatmap_path = paste0(heatmap_path, "/")
  }
  if (!dir.exists(heatmap_path)) {
    dir.create(heatmap_path, recursive = TRUE)
  }
  ## Get the clusters names of the to-test column
  clusters = unique(single.cell.obj@meta.data[, column.to.chisq])
  ## sort the clusters
  clusters = sort(clusters)
  ## Get the meta.data
  meta.data = single.cell.obj@meta.data
  # split the data.frame by the compare.groups and column.to.chisq
  split_df = split(meta.data, list(meta.data[, compare.groups], meta.data[, column.to.chisq]))
  # apply the length function to get counts
  counts = sapply(split_df, nrow)
  # unlist and create a new data.frame
  a = data.frame(do.call(rbind, strsplit(names(counts), "\\.")), a = counts, row.names = NULL)
  # rename the columns
  names(a) = c(compare.groups, column.to.chisq, "a")
  ## Repeat the above steps to for other forms to conduct the chi-square test
  split_df = split(meta.data, list(meta.data[, compare.groups]))
  counts = sapply(split_df, nrow)
  ac = data.frame(do.call(rbind, strsplit(names(counts), "\\.")), ac = counts, row.names = NULL)
  names(ac) = c(compare.groups, "ac")
  split_df = split(meta.data, list(meta.data[, column.to.chisq]))
  counts = sapply(split_df, nrow)
  ab = data.frame(do.call(rbind, strsplit(names(counts), "\\.")), ab = counts, row.names = NULL)
  names(ab) = c(column.to.chisq, "ab")
  d1 = merge(a,ac)
  d1 = merge(d1,ab)
  
  d1["b"] = d1$ab - d1$a
  d1["c"] = d1$ac - d1$a
  d1["d"] = sum(d1$a) - d1$ab - d1$c
  
  ## Conduct the chi-square test
  d1[,c("p.value","Ea","Ec","Eb","Ed")] = t(apply(d1,1,function(x){
    x = as.numeric(x[c('a','c','b','d')])
    k.test = chisq.test(matrix(x, nrow = 2, ncol = 2))
    return(c(k.test$p.value, as.vector(k.test$expected)))
  }))
  d1["Reo"] = d1$a / d1$Ea
  ## Create p.value symbol column
  d1$p.value.symbol = cut(d1$p.value, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c('****', '***', '**', '*', 'N.S.'))
  
  ## Prepare to plot the heatmap:
  require(reshape2)
  require(pheatmap)
  require(ComplexHeatmap)
  dist.heatmap.data = d1[,c(column.to.chisq, compare.groups, "Reo")]
  dist.heatmap = dcast(dist.heatmap.data, dist.heatmap.data[,1] ~ dist.heatmap.data[,2], value.var = "Reo")
  # Now set the row names and remove the column.to.chisq column
  names(dist.heatmap)[1] = column.to.chisq
  rownames(dist.heatmap) = dist.heatmap[,column.to.chisq]
  dist.heatmap = dist.heatmap[, -1]
  ## Round the values for better visualization
  dist.heatmap1 = round(dist.heatmap,digits = 2)
  ## Reorder the heatmap by the cluster order
  dist.heatmap1 = dist.heatmap1[order(factor(rownames(dist.heatmap1), levels = clusters)),]
  
  ## Prepare p-value symbols for heatmap
  dist.heatmap.pval.symbols = d1[, c(column.to.chisq, compare.groups, "p.value.symbol")]
  dist.heatmap.pval.symbols = dcast(dist.heatmap.pval.symbols, dist.heatmap.pval.symbols[,1] ~ dist.heatmap.pval.symbols[,2], value.var = "p.value.symbol")
  names(dist.heatmap.pval.symbols)[1] = column.to.chisq
  rownames(dist.heatmap.pval.symbols) = dist.heatmap.pval.symbols[,column.to.chisq]
  dist.heatmap.pval.symbols = dist.heatmap.pval.symbols[,-1]
  dist.heatmap.pval.symbols = dist.heatmap.pval.symbols[order(factor(rownames(dist.heatmap.pval.symbols), levels = clusters)),]
  
  if (using_complexheatmap) {
    # Create the row annotation for the p-value symbols
    row_anno <- rowAnnotation(p.value = anno_text(dist.heatmap.pval.symbols[,1], 
                                                  just = "center",
                                                  gp = gpar(fontsize = fontsize)))
    # Plot the heatmap
    palette = colorRamp2(c(0, max(dist.heatmap1)/2, max(dist.heatmap1)), c("#FFFFFF",'#FFCC66',"#CC3300"))
    ht = Heatmap(dist.heatmap1,
                 name = "Ro/e",
                 col = palette,
                 show_row_names = TRUE, 
                 show_column_names = TRUE,
                 row_title = column.to.chisq, 
                 column_title = compare.groups, 
                 left_annotation = row_anno, cluster_columns = F, cluster_rows = F, 
                 row_title_gp = gpar(fontsize = 0), 
                 row_names_gp = gpar(fontsize = fontsize), # set fontsize for row names
                 column_names_gp = gpar(fontsize = fontsize), # set fontsize for column names
                 column_title_gp = gpar(fontsize = 0),
                 column_names_rot = 45,
                 width = unit(cellwidth/1.5, "mm") * ncol(dist.heatmap1), # set width of heatmap
                 height = unit(cellheight/2.5, "mm") * nrow(dist.heatmap1),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(round(dist.heatmap1[i, j], 2), x, y, 
                             gp = gpar(fontsize = fontsize))
                 },
                 heatmap_legend_param = list(
                   direction = "vertical", 
                   title_gp = gpar(fontsize = fontsize), # fontsize for legend title
                   labels_gp = gpar(fontsize = fontsize-4),  # fontsize for legend labels
                   legend_height = unit(4, "cm"),
                   grid_width = unit(.8, "cm")
                 )) # set height of heatmap)
    # Save the heatmap
    pdf(file = paste0(heatmap_path, column.to.chisq, "_chisq_test_heatmap.pdf"), width = width, height = height)
    draw(ht, heatmap_legend_side = "right")
    # Add main title
    grid::grid.text(paste0("Chi-square test for ", column.to.chisq, " and ", compare.groups, " in each cluster"), 
                    x = unit(0.5, "npc"), y = unit(0.8, "npc"), just = "center", gp = gpar(fontsize = fontsize))
    dev.off()
  } else {
    if (display_pval_in_heatmap) {
      ## Plot the heatmap
      cell_cluster_dist1 = pheatmap(dist.heatmap1,
                                    cluster_rows = FALSE,
                                    cluster_cols = FALSE,
                                    treeheight_row = 0,
                                    color = heatmap_palette,
                                    display_numbers = dist.heatmap.pval.symbols,
                                    filename = paste0(heatmap_path, column.to.chisq,"_chisq_test_heatmap.pdf"),
                                    cellwidth = cellwidth, cellheight = cellheight, fontsize = fontsize, width = width, height = height,
                                    border_color = NA,
                                    angle_col='45',
                                    main = paste0("Chi-square test for ", column.to.chisq, " and ", compare.groups, " in each cluster\n"),
                                    angle_row='45')
    } else {
      p.anno = dist.heatmap.pval.symbols[,1, drop = F]
      names(p.anno) = "p.value"
      cell_cluster_dist1 = pheatmap(dist.heatmap1,
                                    cluster_rows = FALSE,
                                    cluster_cols = FALSE,
                                    treeheight_row = 0,
                                    color = heatmap_palette,
                                    display_numbers = dist.heatmap.pval.symbols,
                                    filename = paste0(heatmap_path, column.to.chisq,"_chisq_test_heatmap.pdf"),
                                    cellwidth = cellwidth, cellheight = cellheight, fontsize = fontsize, width = width, height = height,
                                    border_color = NA,
                                    angle_col='45',
                                    main = paste0("Chi-square test for ", column.to.chisq, " and ", compare.groups, " in each cluster\n"),
                                    angle_row='45', annotation_row = p.anno)
    }
  }
}