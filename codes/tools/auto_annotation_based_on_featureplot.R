auto_annotation_based_on_featureplot = function(genesets, scRNA, violin = F, combine_fea_vio = F,
                                                outdir = "cell_identify/single_cell_marker/"){
  scmarkeroutput = function(scRNA, features, cellname, violin = F, combine_fea_vio = F,
                            outdir){
    set.seed(123)
    if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
    tmp = FeaturePlot(scRNA, features = features, reduction = "umap", label = T, pt.size = .1)
    # The following script decides whether to generate a violin plot：
    if (violin){
      vln = VlnPlot(scRNA, features = features, pt.size = 0.2, group.by = "seurat_clusters") +
        scale_fill_discrete(name = "Cluster")
      if (combine_fea_vio){
        v.f.plot = tmp|vln
        ggsave(paste0(outdir,cellname,"_combine_plot.pdf"), plot = v.f.plot, width = 12,
               height = 8)
      } else {
        ggsave(paste0(outdir,cellname,"_vlnplot.pdf"), plot = vln, width = 12, height = 8)
      }
    } else {
      ggsave(paste0(outdir,cellname,"_feaplot.pdf"), plot = tmp, width = 12, height = 8)
    }
  }
  ## Check if the genesets contains "", or NA and remove them
  if (is.null(genesets)) {
    stop("genesets is NULL")
  } else if (nrow(genesets) == 0) {
    stop("genesets is empty")
  } else if ("" %in% names(genesets)) {
    genesets = genesets[,-which(names(genesets) == "")]
    genesets = genesets %>% as.list()
    if ("" %in% genesets | any(is.na(genesets))) {
      genesets = lapply(genesets, function(x) x[!is.na(x) & x != ""])
    } else {
      genesets = lapply(genesets, function(x) x[x != ""])
    }
  }
  ## Generate the marker output
  for (cells in names(genesets)){
    cat("Plotting", cells, "FeaturePlots...\n")
    genes = genesets[[cells]]
    scmarkeroutput(scRNA, features = genes, cellname = cells, violin = violin,
                   combine_fea_vio = combine_fea_vio, outdir = outdir)
  }
}
