dimplot_inhouse_R <- function(obj, plot_item, plot_order = NULL, prefix = "BGI_output/default", point_size = 0.55,
                              plot_region = F, bg_color = "#BDBDBD"){
  require(tidyverse)
  require(glue)
  write.table(obj@images$slice1@coordinates, 'df_pos_tmp.txt', quote=F, sep='\t')
  write.table(obj@meta.data, 'df_meta_tmp.txt', quote=F, sep='\t')
  tmp = unique(obj@meta.data[,plot_item]) %>% .[order(.)] %>% str_c(., collapse=' ')
  python_path <- system("which python", intern = T)
  script_path <- 'BGI_ST_Dimplot_Featureplot_CJY_for_annotation.py'
  cols = length(unique(obj@meta.data[,plot_item]))
  command <- glue("{python_path}  {script_path} -i {plot_item} --plot_type DimPlot -o {prefix} -r {plot_region} -b '{bg_color}' --plot_order '{tmp}' --length {cols} --point_size {point_size}")
  system(command)
  system('rm df_pos_tmp.txt df_meta_tmp.txt')
}
