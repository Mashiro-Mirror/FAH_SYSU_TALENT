if(T){
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggsci)
  library(png)
  library(glue)
  library(RColorBrewer)
  library(tibble)
  library(tidyr)
  library(harmony)
  library(cowplot)
  library(scRepertoire)
  library(circlize)
  library(scales)
  library(RColorBrewer)
  library(gplots)
  library(batchelor)
  library(Matrix)
  library(ggpubr)
  library(fgsea)
  library(reshape2)
  library(viridis)
  library(ggrepel)
}

featureplot_inhouse <- function(obj, plot_item, prefix){
  write.table(obj@images$slice1@coordinates, 'df_pos_tmp.txt', quote=F, sep='\t')
  write.table(obj@meta.data, 'df_meta_tmp.txt', quote=F, sep='\t')
  command <- glue("/../../../anaconda3/bin/python3.8 {script_path} -i {plot_item} --plot_type FeaturePlot -o {prefix}")
  system(command)
  system('rm df_pos_tmp.txt df_meta_tmp.txt')
}


###############################################################################
#'                   Manuscipt: Figure S2                                    '#
###############################################################################

Sample_list <- c("pt03","pt05","pt10","pt17")
script_path <- './Visualizaton.py'
for (i in Sample_list) {
  obj <- readRDS(glue("./ST{i}.rds"))
  featureplot_inhouse(obj, 'nFeature_Spatial', glue("./{i}_nFeature"))
  featureplot_inhouse(obj, 'nCount_Spatial', glue("./{i}_nCount"))
}















