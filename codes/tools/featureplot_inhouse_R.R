featureplot_inhouse_R <- function(obj, feature, prefix = "BGI_output/default", name = ""){
  require(tidyverse)
  require(stringr)
  require(tidyverse)
  require(Seurat)
  require(reticulate)
  require(glue)
  if(!dir.exists(prefix)){
    dir.create(prefix, recursive = T)
  }
  # Write coordinates to a temporary file
  write.table(obj@images$slice1@coordinates, glue('{prefix}/df_pos_tmp.txt'), quote=F, sep='\t')
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
  write.table(meta.data, glue('{prefix}/df_meta_tmp.txt'), quote=F, sep='\t')
  rm(exp, fea, meta.data)
  gc()
  
  python_path = system("which python", intern = T)
  script_path <- 'tools/plot_cluster.py'
  
  # Prepare the features string
  # feature = paste(feature, collapse = ",")
  command <- glue("{python_path} {script_path} --plot_type FeaturePlot -i {feature} -o {prefix} -n {name}")
  system(command)
  # colorbar
  command <- glue("{python_path} {script_path} -i {feature} --plot_type Colorbar -o {prefix} -n {name}")
  system(command)
  system(glue('rm {prefix}/df_pos_tmp.txt {prefix}/df_meta_tmp.txt'))
}