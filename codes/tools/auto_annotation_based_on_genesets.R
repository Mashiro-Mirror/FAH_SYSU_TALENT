auto_annotation_based_on_genesets = function(genesets,scRNA){
  require(AUCell)
  require(Seurat)
  require(Matrix)
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
  ## AUCell Method
  ## Get the expression matrix and make it as a sparse matrix
  exprMatrix = as(GetAssayData(scRNA, slot = "counts"), "sparseMatrix")
  ## Perform AUCell analysis
  for (geneset in names(genesets)){
    genes = genesets[[geneset]]
    cell_rank = AUCell_buildRankings(exprMatrix, nCores = 12, plotStats = T)
    cell.auc = AUCell_calcAUC(genes, cell_rank)
    auc_value = as.numeric(getAUC(cell.auc)[1,])
    ## Add the auc value to the new column of metadata in the end of the Seurat meta.data
    scRNA@meta.data[,(ncol(scRNA@meta.data)+1)] = auc_value
    names(scRNA@meta.data)[ncol(scRNA@meta.data)] = paste0(geneset,"_AUCell")
  } 
  
  ## AddModuleScore Method
  scRNA = AddModuleScore(scRNA, features = genesets)
  names(scRNA@meta.data)[c((ncol(scRNA@meta.data)-length(genesets)+1):ncol(scRNA@meta.data))] = names(genesets)
  return(scRNA)
}