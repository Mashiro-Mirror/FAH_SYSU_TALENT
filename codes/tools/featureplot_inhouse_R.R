featureplot_inhouse_R <- function(obj, feature, prefix = "BGI_output/default", point_size = 0.55,
                                  num_cols = 4){
  require(tidyverse)
  require(stringr)
  require(tidyverse)
  require(Seurat)
  require(reticulate)
  require(glue)
  # Write coordinates to a temporary file
  write.table(obj@images$slice1@coordinates, 'df_pos_tmp.txt', quote=F, sep='\t')
  # Get the assay data and transform it into a data frame
  exp = GetAssayData(obj, slot = "counts", assay = "Spatial") %>% as_matrix()
  if (length(feature) <= 1) {
    fea = as.data.frame(exp[feature,])
    names(fea) = feature
  } else {
    fea = as.data.frame(t(exp[feature,]))
  }
  meta.data = cbind(obj@meta.data, fea)
  # Write metadata to a temporary file
  write.table(meta.data, 'df_meta_tmp.txt', quote=F, sep='\t')
  rm(exp, fea, meta.data)
  gc()
  
  python_path = system("which python", intern = T)
  script_path <- 'BGI_ST_Dimplot_Featureplot_CJY_for_annotation.py'
  
  # Prepare the features string
  feature = paste(feature, collapse = ",")
  command <- glue("{python_path} {script_path} -f {feature} --plot_type FeaturePlot -o {prefix} --point_size {point_size} --n_columns {num_cols}")
  system(command)
  system('rm df_pos_tmp.txt df_meta_tmp.txt')
}