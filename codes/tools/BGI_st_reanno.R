st_anno_score <- function(obj_merged, gl) {
  obj_merged <- NormalizeData(obj_merged)
  obj_merged <- AddModuleScore(obj_merged, gl)
  mtx_score <- obj_merged@meta.data %>%
    dplyr::select(all_of(setNames(str_glue("Cluster{seq_along(gl)}"), names(gl)))) 
  return(mtx_score)
}
#####对bin进行分组，每组bin不超过65000
st_anno_split_group <- function(mtx_score, group_size = 60000, seed = 123) {
  n_group <- ceiling(nrow(mtx_score) / group_size)
  set.seed(seed)
  dat_group <- tibble(
    id = rownames(mtx_score), 
    rand = runif(nrow(mtx_score)), 
    group = cut_number(rand, n_group) %>% as.numeric() 
  )
  split(dat_group$id, dat_group$group)
}
#####限速步骤，对bin进行层次聚类，对每个聚类取均值，分配属于哪一个cluster
st_anno_hclust_bin <- function(mtx_score_split, mtx_score_scale) {
  # Hierarchical Clustering
  tree_bin <- pheatmap:::cluster_mat(as.matrix(mtx_score_scale), distance = "euclidean", method = "complete")
  k_init <- min(200, nrow(mtx_score_scale))
  bin_cluster <- cutree(tree_bin, k = k_init)
  
  # Annotation with Mean Score of Cluster
  df_anno <- data.frame(bin = names(bin_cluster), cluster = bin_cluster)
  df_cluster_score <- bind_cols(df_anno, mtx_score_scale[rownames(df_anno), ,drop = F]) %>% 
    pivot_longer(all_of(colnames(mtx_score)), names_to = "celltype", values_to = "score") %>% 
    # compute the mean of each cluster
    group_by(cluster, celltype) %>% 
    dplyr::summarise(score = sum(score) / n_distinct(bin), .groups = "drop") %>% 
    # determine cluster by maximum value
    group_by(cluster) %>% 
    dplyr::filter(score == max(score)) %>% 
    # bins with multiple maximum scores were "unknown"
    group_by(cluster) %>% 
    dplyr::mutate(n = n()) %>% 
    ungroup() 
  df_bin_score <- df_cluster_score %>% 
    left_join(df_anno, by = "cluster", multiple = "all")
  return(df_bin_score)
}
#####设定unknown阈值，给每个bin确定一类细胞
st_anno_score_cut <- function(bin_score, mtx_score_scale, prob = NULL, value = NULL) {
  if (is.null(prob) + is.null(value) == 1) {
    if (!is.null(prob)) {
      cutoff_score <- bin_score %>% distinct(cluster, score) %>% pull(score) %>% quantile(prob)
    }
    if (!is.null(value)) {
      cutoff_score <- value
    }
  }
  bin_subtype <- bin_score %>% dplyr::filter(n == 1, score > cutoff_score) 
  bin_unknown <- bin_score %>% filter(!bin %in% bin_subtype$bin) %>% mutate(celltype = "unknown")
  df_subtype <- bind_rows(bin_subtype, bin_unknown) %>% dplyr::select(bin, celltype, score)
  
  dat_p <- bind_cols(df_subtype, mtx_score_scale[df_subtype$bin, ,drop = F]) %>% 
    pivot_longer(all_of(colnames(mtx_score_scale)))
  return(dat_p)
}

st_anno_plot_score <- function(bin_score_celltype, gl) {
  lvl_celltype <- c(names(gl), "unknown")
  lvl_bin <- bin_score_celltype %>% 
    mutate(across(celltype, ~ factor(.x, levels = lvl_celltype))) %>% 
    group_split(celltype) %>% 
    map(~ {
      .x %>% 
        filter(celltype == name | celltype == "unknown") %>% 
        arrange(desc(value)) %>% 
        pull(bin) %>% 
        unique()
    }) %>% 
    unlist()
  q <- quantile(bin_score_celltype$value, probs = c(0.05, 0.95))
  bin_score_celltype %>% 
    mutate(
      across(bin, ~ factor(.x, levels = lvl_bin)), 
      across(c(name, celltype), ~ factor(.x, levels = lvl_celltype)), 
      across(value, ~ ifelse(.x > q["95%"], q["95%"], .x)), 
      across(value, ~ ifelse(.x < q["5%"], q["5%"], .x))
    ) %>% 
    ggplot(aes(x = bin, y = name, fill = value)) +
    geom_tile() +
    facet_grid(~ celltype, scales = "free", space = "free") +
    labs(x = NULL, y = NULL, fill = NULL) +
    scale_fill_gradientn(
      colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(256)
    ) +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()
    )
}