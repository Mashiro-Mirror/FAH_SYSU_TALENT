dimplot_inhouse_R <- function(obj, plot_item, tmp_color = ggsci::pal_igv()(50), plot_order = NULL, prefix = "BGI_output/default"){
  require(tidyverse)
  require(glue)
  if(!dir.exists(prefix)){
    dir.create(prefix, recursive = T)
  }
  write.table(obj@images$slice1@coordinates, glue('{prefix}/df_pos_tmp.txt'), quote=F, sep='\t')
  write.table(obj@meta.data, glue('{prefix}/df_meta_tmp.txt'), quote=F, sep='\t')
  tmp = unique(obj@meta.data[,plot_item]) %>% .[order(.)]
  if ("Others" %in% tmp) {
    tmp = tmp[-which(tmp == "Others")]
  }
  tmp = tmp %>% str_c(., collapse=' ')
  
  python_path <- system("which python", intern = T)
  script_path <- 'tools/plot_cluster.py'
  cols = length(unique(obj@meta.data[,plot_item]))
  tmp_color <- str_c(tmp_color, collapse=' ')
  command <- glue("{python_path}  {script_path} -i {plot_item} --plot_type DimPlot -o {prefix} --plot_order '{tmp}' --plot_color '{tmp_color}' -n {plot_item}")
  system(command)
  system(glue('rm {prefix}/df_pos_tmp.txt {prefix}/df_meta_tmp.txt'))
}
