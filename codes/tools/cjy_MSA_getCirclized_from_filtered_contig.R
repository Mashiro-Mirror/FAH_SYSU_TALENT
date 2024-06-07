# Compute total within-cluster sum of square
tot.withinss <- function(cluster, dist_matrix) {
  cluster_ids <- levels(as.factor(cluster))
  withinss <- numeric(length(cluster_ids))
  
  for (i in seq_along(cluster_ids)) {
    cluster_points <- which(cluster == cluster_ids[i])
    withinss[i] <- sum(dist_matrix[cluster_points, cluster_points])
  }
  
  sum(withinss)
}

cjy_MSA_from_filtered_contig = function(merged_contigs,number_range_of_clonotype=500) {
  require(reshape2)
  require(msa)
  require(circlize)
  require(ggpubr)
  
  ## Perform MSA and re cluster the cdrs to assign them to different clonotypes
  # Create a AAStringSet object
  sequences = Biostrings::AAStringSet(merged_contigs$cdrs)
  # Perform multiple sequence alignment
  alignment <- msa(sequences)
  # Convert alignment to character matrix
  alignment_matrix <- as.character(alignment)
  
  # Compute pairwise distance matrix
  dist_matrix <- stringdist::stringdistmatrix(alignment_matrix, alignment_matrix)
  # Perform hierarchical clustering
  cluster <- hclust(as.dist(dist_matrix))
  # Compute wss for k = 2 to k = number_range_of_clonotype (default = 500)
  k.values <- 2:number_range_of_clonotype
  wss.values <- sapply(k.values, function(k) tot.withinss(cutree(cluster, k), dist_matrix))
  # Compute the rate of decrease of WSS
  rate_of_decrease <- diff(wss.values)
  # Determine the 'elbow point'
  optimal_clusters <- which.max(rate_of_decrease) + 1
  # Cut the hierarchical clustering tree with the optimal number of clusters
  clusters <- cutree(cluster, optimal_clusters)
  
  # Create a data frame with sequences and their corresponding cluster ID
  result <- data.frame(cdrs = merged_contigs$cdrs, 
                       new_Clonotype_ID = paste0("new_clonotype",clusters))
  ## Merge with the merged_contigs
  new_merged_contigs = cbind(merged_contigs, result[,2])
  names(new_merged_contigs)[ncol(new_merged_contigs)] = "new_Clonotype_ID"
  return(new_merged_contigs)
}

cjy_Identical_cdr3nt_from_filtered_contig = function(merged_contigs) {
  ## From the cdr3nt column, we can get the identical cdr3nt and assign them to a new clonotype
  ## First, we need to get the identical cdr3nt
  identical_cdr3nt = merged_contigs %>% group_by(cdr3_nt) %>% summarise(Freq = n())
  identical_cdr3nt$new_Clonotype_ID = paste0("new_clonotype",1:nrow(identical_cdr3nt))
  ## Merge with the merged_contigs
  new_merged_contigs = merge(merged_contigs, identical_cdr3nt, by = "cdr3_nt")
  return(new_merged_contigs)
}

cjy_percent_Identical_cdr3nt_from_filtered_contig = function(merged_contigs, percent_identical = 75) {
  require(stringdist)
  require(dplyr)
  ## Convert percent_identical to a decimal
  percent_identical_decimal = percent_identical / 100
  ## Calculate the Jaccard distance between all pairs of cdr3_nt sequences
  cdr3_nt = merged_contigs$cdr3_nt
  dist_matrix = stringdist::stringdistmatrix(cdr3_nt, cdr3_nt, method = "jaccard")
  ## Identify pairs of cdr3_nt sequences that are at least percent_identical_decimal identical
  ## (Jaccard distance <= 1 - percent_identical_decimal, because it's a measure of dissimilarity)
  identical_pairs = which(dist_matrix <= 1 - percent_identical_decimal, arr.ind = TRUE)
  ## Remove pairs where the cdr3_nt sequences are different lengths
  identical_pairs = identical_pairs[nchar(cdr3_nt[identical_pairs[, 1]]) == nchar(cdr3_nt[identical_pairs[, 2]]), ]
  ## Create a list of cdr3_nt groups, where each group contains cdr3_nt sequences that are at least percent_identical_decimal identical
  cdr3_nt_groups = list()
  assigned_cdr3_nt = character(0)
  for (pair in seq_len(nrow(identical_pairs))) {
    if (!(identical_pairs[pair, 1] %in% assigned_cdr3_nt)) {
      cdr3_nt_groups[[length(cdr3_nt_groups) + 1]] = identical_pairs[pair, 1]
      assigned_cdr3_nt = c(assigned_cdr3_nt, identical_pairs[pair, 1])
    }
  }
  ## Create a new clonotype ID for each group
  new_Clonotype_ID = unlist(lapply(seq_along(cdr3_nt_groups), function(i) rep(paste0("new_clonotype", i), length(cdr3_nt_groups[[i]]))))
  ## Create a data frame mapping each cdr3_nt to its new clonotype ID
  new_clonotype_df = data.frame(cdr3_nt = cdr3_nt[unlist(cdr3_nt_groups)], new_Clonotype_ID = new_Clonotype_ID)
  ## Remove duplicates in new_clonotype_df
  new_clonotype_df = new_clonotype_df[!duplicated(new_clonotype_df$cdr3_nt), ]
  ## Merge with the merged_contigs
  new_merged_contigs = merge(merged_contigs, new_clonotype_df, by = "cdr3_nt")
  return(new_merged_contigs)
}


cjy_getCirclized_after_cjy_MSA = function(new_merged_contigs, clonocall = "new_Clonotype_ID", 
                                          group.by = "orig.ident", plot_dir, sample, 
                                          n_samples = 400,
                                          plot.width = 8, plot.height = 8,
                                          chord_transparency = .5,
                                          outer.grid.space = .5,
                                          link.not.display.threshold = 0.1,
                                          niceFacing = T, if_from_to = F,
                                          big.gap = .5, small.gap = .2) {
  require(tidyverse)
  require(RColorBrewer)
  require(circlize)
  ## judge if the new_merged_contigs is empty
  if (nrow(new_merged_contigs) == 0) {
    cat("The new_merged_contigs is empty, please check the input data!\n")
    return()
  }
  ## Make sure the cloneCall is in the right scRepertoire format
  cloneCall <- clonocall
  ## If group.by is NULL, then use the ident column
  test <- new_merged_contigs[, c(cloneCall, group.by)]
  ## Remove NA clonotype calls
  test <- test[!is.na(test[, cloneCall]), ]
  test$clone_source = paste0(test[, group.by], "_", test[, cloneCall])
  ## Create a blank matrix to store the custom data
  unique.new.clonotypes <- unique(test[, "clone_source"])
  if (length(unique.new.clonotypes) == 1) {
    cat("Only one unique new clonotype, please check the input data!\n")
    return()
  }
  cloneCallData <- matrix(0, nrow = length(unique.new.clonotypes), 
                          ncol = length(unique.new.clonotypes))
  ## Assign rownames and colnames to the matrix
  rownames(cloneCallData) <- unique.new.clonotypes
  colnames(cloneCallData) <- unique.new.clonotypes
  ## Fill in the matrix with the following algorithm
  for (unique.new.clonotype in unique(test[, cloneCall])) {
    ## Get the current clone and the corresponding orig.ident
    current_clone <- test[test[, cloneCall] == unique.new.clonotype, ]
    ## Convert the current_clone to table
    clone_table <- table(current_clone[,cloneCall], current_clone[, group.by]) %>% as.data.frame()
    ## Judge if the clone_table has more than 1 row
    if (nrow(clone_table) > 1) {
      cat(unique.new.clonotype,"\n")
      ## list all the combinations of the orig.ident (Var2)
      combinations <- combn(unique(clone_table[,2]), 2)
      ## increment the corresponding cell in the matrix by 1 for each combination
      for (i in seq_len(ncol(combinations))) {
        comb <- combinations[, i]
        freq <- clone_table[clone_table[,2] %in% comb, "Freq"]
        ## Here we need to judge if the combination is already in the rownames and !! columns names of matrix
        if (paste0(comb[1], "_", unique.new.clonotype) %in% rownames(cloneCallData) && paste0(comb[2], "_", unique.new.clonotype) %in% rownames(cloneCallData)) {
          if (if_from_to) {
            cloneCallData[paste0(comb[1], "_", unique.new.clonotype),
                          paste0(comb[2], "_", unique.new.clonotype)] <- cloneCallData[paste0(comb[1], "_", unique.new.clonotype),
                                                                                       paste0(comb[2], "_", unique.new.clonotype)] + freq[2]
            cloneCallData[paste0(comb[2], "_", unique.new.clonotype),
                          paste0(comb[1], "_", unique.new.clonotype)] <- cloneCallData[paste0(comb[2], "_", unique.new.clonotype),
                                                                                       paste0(comb[1], "_", unique.new.clonotype)] + freq[1]
          } else {
            cloneCallData[paste0(comb[1], "_", unique.new.clonotype),
                          paste0(comb[2], "_", unique.new.clonotype)] <- cloneCallData[paste0(comb[1], "_", unique.new.clonotype),
                                                                                       paste0(comb[2], "_", unique.new.clonotype)] + sum(freq)
          }
        } else {
          next
        }
      }
    } else {
      ## The similar judgement as above
      if (paste0(clone_table[,2], "_", unique.new.clonotype) %in% rownames(cloneCallData)) {
        freq <- clone_table[1, "Freq"]
        cloneCallData[paste0(clone_table[,2], "_", unique.new.clonotype),
                      paste0(clone_table[,2], "_", unique.new.clonotype)] <- freq
      } else {
        next
      }
    }
  }
  # Add a small non-zero value to rows with all zero values for chordDiagram plotting 
  # cloneCallData[cloneCallData == 0] <- 1e-3
  ## To avoid self-link, set the diagonal to 0
  diag(cloneCallData) <- 0
  
  if (nrow(cloneCallData) > n_samples) {
    # Get the total number of rows/columns in the matrix
    n_total <- nrow(cloneCallData)
    # Generate a random sample of row/column indices, allowing for replacement
    set.seed(123)
    sample_indices <- sample(n_total, n_samples, replace = TRUE)
    # Subset the matrix using the sample indices
    cloneCallData <- cloneCallData[sample_indices, sample_indices]
  } else {
    cloneCallData <- cloneCallData
  }
  
  ## Finish the chordDiagram plot
  nm = unique(unlist(dimnames(cloneCallData)))
  group = structure(gsub("_new_clonotype[0-9]*$","", nm), names = nm)
  ## Assign group variable to the group argument and add grid.col for plotting
  grid.freq = as.data.frame(table(group))
  grid.col = c()
  for (i in 1:nrow(grid.freq)) {
    g.col = structure(rep("#FFFFFF", as.numeric(grid.freq[i, "Freq"])), 
                      names = names(group)[group==as.character(grid.freq[i,1])])
    grid.col = c(grid.col, g.col)
  }
  palette_function <- colorRampPalette(c("red"))
  col_mat = palette_function(length(cloneCallData))
  dim(col_mat) = dim(cloneCallData)
  # Add transparency
  alpha_value = chord_transparency  # set the transparency level between 0 (completely transparent) and 1 (completely opaque)
  col_mat <- apply(col_mat, MARGIN = c(1,2), FUN = function(color) grDevices::adjustcolor(color, alpha.f = alpha_value))
  col_mat[cloneCallData < ((max(cloneCallData) - min(cloneCallData)) * link.not.display.threshold) | cloneCallData < 1] <- "#00000000"
  ## make sure the plot_dir exists:
  if(!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  pdf(paste0(plot_dir,sample,"_chordDiagram.pdf"), width = plot.width, height = plot.height)
  chordDiagram(cloneCallData, group = group, grid.col = grid.col, 
               annotationTrack = c("grid"), small.gap = small.gap,
               big.gap = big.gap,
               col = col_mat, transparency = chord_transparency,
               grid.border = "#000000",
               preAllocateTracks = list(track.height = mm_h(6),
                                        track.margin = c(mm_h(outer.grid.space), 0)))
  # ## Here is a difficult part, we need to add the outer track for the grid.
  ## But the number of the grid, to be more accurate, the types of the grid is uncertain.
  ## So we need to judge the types of the grid and add the outer track for each type.
  ## First, we need to get the types of the grid
  grid.types = unique(grid.freq[,1]) %>% as.character()
  grid.pals = brewer.pal(length(grid.types), "Set1")
  pals = 1
  for (unique.type in grid.types){
    print(unique.type)
    grid.names = rownames(cloneCallData)[grep(paste0(unique.type,"_"), rownames(cloneCallData))]
    ## Judge if the length of grid.names is larger than a threshold:
    highlight.sector(grid.names, track.index = 1, col = grid.pals[pals], 
                     text = unique.type, cex = 1, text.col = "white", 
                     niceFacing = niceFacing)
    pals = pals + 1
  }
  circos.clear()
  dev.off()
  ## Reset the value for output
  cloneCallData[cloneCallData == 0.001] = 0
  openxlsx::write.xlsx(as.data.frame(cloneCallData), file = paste0(plot_dir,sample,".xlsx"), rowNames = T)
}

cjy_getCirclized_after_cjy_MSA_from_data_frame = function(clonotype.df, plot_dir, sample,
                                                          plot.width = 8, plot.height = 8,
                                                          chord_transparency = .5,
                                                          outer.grid.space = .5,
                                                          link.not.display.threshold = 0.1,
                                                          niceFacing = T,big.gap = .5, 
                                                          small.gap = .2, seed = 123) {
  ## Make sure that the first 4 columns are: from, to, value1, value2
  if (ncol(clonotype.df) < 4) {
    cat("The input data frame should have at least 4 columns: from, to, value1, value2\n")
    return()
  } else {
    if (!("from" %in% colnames(clonotype.df) && "to" %in% colnames(clonotype.df) && 
          "value1" %in% colnames(clonotype.df) && "value2" %in% colnames(clonotype.df))) {
      cat("The input data frame should have at least 4 columns: from, to, value1, value2\n")
      return()
    }
  }
  ## Finish the chordDiagram plot
  nm = unique(c(clonotype.df$from, clonotype.df$to))
  group = structure(gsub("_new_clonotype[0-9]*$","", nm), names = nm)
  ## Assign group variable to the group argument and add grid.col for plotting
  grid.freq = as.data.frame(table(group))
  grid.col = c()
  for (i in 1:nrow(grid.freq)) {
    g.col = structure(rep("#FFFFFF", as.numeric(grid.freq[i, "Freq"])), 
                      names = names(group)[group==as.character(grid.freq[i,1])])
    grid.col = c(grid.col, g.col)
  }
  palette_function <- colorRampPalette(c("red"))
  set.seed(seed)
  col_mat = rand_color(nrow(clonotype.df))
  # Add transparency
  alpha_value = chord_transparency  # set the transparency level between 0 (completely transparent) and 1 (completely opaque)
  col_mat = sapply(col_mat, function(color) grDevices::adjustcolor(color, alpha.f = alpha_value))
  ## make sure the plot_dir exists:
  if(!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  pdf(paste0(plot_dir,sample,"_chordDiagram_data_frame.pdf"), width = plot.width, height = plot.height)
  chordDiagram(clonotype.df, grid.col = grid.col, 
               annotationTrack = c("grid"), small.gap = small.gap,
               big.gap = big.gap,
               col = col_mat, transparency = chord_transparency,
               grid.border = "#000000",
               preAllocateTracks = list(track.height = mm_h(6),
                                        track.margin = c(mm_h(outer.grid.space), 0)))
  # ## Here is a difficult part, we need to add the outer track for the grid.
  ## But the number of the grid, to be more accurate, the types of the grid is uncertain.
  ## So we need to judge the types of the grid and add the outer track for each type.
  ## First, we need to get the types of the grid
  grid.types = unique(grid.freq[,1]) %>% as.character()
  grid.pals = brewer.pal(length(grid.types), "Set1")
  pals = 1
  for (unique.type in grid.types){
    print(unique.type)
    grid.names = clonotype.df[, grep(paste0(unique.type,"_"), clonotype.df[1,,drop = TRUE])]
    ## Judge if the length of grid.names is larger than a threshold:
    highlight.sector(grid.names, track.index = 1, col = grid.pals[pals], 
                     text = unique.type, cex = 1, text.col = "white", 
                     niceFacing = niceFacing)
    pals = pals + 1
  }
  circos.clear()
  dev.off()
  openxlsx::write.xlsx(clonotype.df, file = paste0(plot_dir,sample,"_data_frame.xlsx"), rowNames = T)
}
