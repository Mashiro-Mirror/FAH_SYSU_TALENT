# Define a function to extract and concatenate the parts of the input vector
extract_parts <- function(input_vec, symbol = "-", part.to.choose) {
  require(Seurat)
  require(tidyverse)
  require(openxlsx)
  require(ggpubr)
  require(glue)
  require(jsonlite)
  # Extract parts of each string
  result <- t(sapply(input_vec, function(s) {
    split_string <- strsplit(s, symbol)
    c(part1 = split_string[[1]][1], part2 = split_string[[1]][2])
  }))
  
  parts = result[, part.to.choose]
  
  return(parts)
}

interactionIndexGeneration = function(stRNA = stRNA, max_distance = 1, sender_Cells, receiver_Cells, Regions = "Tu",
                                      outputDirectory, slot = "counts", ligand, receptor, bg_color = "#BDBDBD"){
  stRNA = stRNA
  stRNA$only_bins = gsub(paste0("_",unique(stRNA$orig.ident)),"", stRNA$bins)
  meta.data = stRNA@meta.data
  after.meta = stRNA@meta.data
  after.meta$Interaction_Pair = "Others"
  ##get sender cells coordinate
  interaction.set = subset(stRNA, new_cellsubtype == sender_Cells)
  interaction.meta = interaction.set@meta.data
  interaction.bins = interaction.meta$only_bins
  interaction.coords = interaction.meta[,c("only_bins","row","col","new_cellsubtype")]
  after.meta$Interaction_Pair[after.meta$only_bins %in% interaction.bins] = sender_Cells
  ##get surrounding receiver bins:
  surround.coords = get_all_surrounding_coords_withBins(interaction.coords,max_distance = max_distance)
  surround.meta = merge(after.meta, surround.coords, by = c("row","col"))
  surround.bins = surround.meta$only_bins[!surround.meta$only_bins %in% interaction.bins & surround.meta$new_cellsubtype == receiver_Cells & surround.meta$Bin_Region == Regions]
  after.meta$Interaction_Pair[after.meta$bins %in% surround.bins] = paste0(receiver_Cells)
  stRNA@meta.data = after.meta
  ##Plotting BGI Dimplot
  plot_item = "Interaction_Pair"
  prefix = outputDirectory
  ## check the column exist in meta.data
  if (plot_item %in% names(stRNA@meta.data)) {
    dimplot_inhouse_R(obj = stRNA, plot_item = plot_item, prefix = prefix, bg_color = bg_color)
  } else {
    tryCatch({
      print(paste0(id, " ",plot_item," not in the meta.data!"))
    }, error = function(e) {
      print(paste0("Error occurred in iteration for ", id, ": ", e$message))
    })
  }
  ## Get Sender-Receiver Pair Data.frame
  if (length(interaction.bins) > 0 & length(surround.bins) > 0){
    only.meta = surround.meta[!surround.meta$only_bins %in% interaction.bins & surround.meta$new_cellsubtype == receiver_Cells & surround.meta$Bin_Region == Regions,]
    only.meta$SR_Bins = paste0(only.meta$senderBins, "-", only.meta$only_bins)
    only.meta$SR_Type = paste0(only.meta$senderType, "-", only.meta$new_cellsubtype)
    ## Get expression data
    expr.bins = c(only.meta$only_bins, only.meta$senderBins)
    sub.st = subset(stRNA, only_bins %in% expr.bins)
    expr = GetAssayData(sub.st, slot = slot, assay = "Spatial") %>% as_matrix()
    colnames(expr) = gsub(paste0("_",unique(stRNA$orig.ident)),"", colnames(expr))
    ## Generate the LR-interaction potential matrix:
    lr.df = data.frame(row.names = only.meta$SR_Bins,
                       lr_pair = rep("", length(only.meta$SR_Bins)))
    colnames(lr.df)[1] = paste0(ligand, "-",receptor)
    for (pair_num in 1:nrow(only.meta)){
      split_string = strsplit(only.meta$SR_Bins[pair_num],"-")[[1]]
      sender_bins <- split_string[1]
      receiver_bins <- split_string[2]
      ligand.exp = expr[ligand,sender_bins]
      receptor.exp = expr[receptor,receiver_bins]
      interaction.potential = atan(ligand.exp * receptor.exp) / (pi / 2)
      # interaction.potential = ligand.exp * receptor.exp
      lr.df[pair_num,1] = interaction.potential
    }
    lr.df[,1] = as.numeric(lr.df[,1])
    # lr.df[,1] = (lr.df[,1] - min(lr.df[,1])) / (max(lr.df[,1]) - min(lr.df[,1]))
    return(lr.df)
  } else {
    tryCatch({
      cat(unique(stRNA$orig.ident), "doesn't have either current sender or receiver.\n")
    }, error = function(e) {
      cat("Error occurred in ", unique(stRNA$orig.ident), ": ", e$message,"\n")
    })
  }
  
}

dist2list <- function(dist){
  if(!class(dist) == "dist"){
    stop("the input data must be a dist object.")
  }
  dat <- as.data.frame(as.matrix(dist))
  if(is.null(names(dat))){
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames,rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)
  return(res)
}

interactionDistanceGeneration = function(stRNA = stRNA, sender_Cells, receiver_Cells, Regions = "Tu",
                                         ligand, receptor, bg_color = "#BDBDBD"){
  stRNA = stRNA
  stRNA$only_bins = gsub(paste0("_",unique(stRNA$orig.ident)),"", stRNA$bins)
  id = unique(stRNA$orig.ident)
  ##get ligand-positive sender cells coordinate
  sender.set = NULL
  receiver.set = NULL
  tryCatch({
    sender.set = subset(stRNA, new_cellsubtype == sender_Cells & Bin_Region == Regions & FetchData(object = stRNA, vars = ligand, slot = "counts") > 0)
  }, error = function(e){
    cat("Error encountered:", e$message, "\n")
    return(NULL)
  })
  if (!is.null(sender.set)) {
    sender.df = sender.set@meta.data[,c("only_bins","row","col", "new_cellsubtype")]
    sender.bins = rownames(sender.df)
  } else {
    return(NULL)
  }
  
  ##get receptor-positive receiver cells coordinate
  tryCatch({
    receiver.set = subset(stRNA, new_cellsubtype == receiver_Cells & Bin_Region == Regions & FetchData(object = stRNA, vars = receptor, slot = "counts") > 0)
  }, error = function(e){
    cat("Error encountered:", e$message, "\n")
    return(NULL)
  })
  if (!is.null(receiver.set)) {
    receiver.df = receiver.set@meta.data[,c("only_bins","row","col", "new_cellsubtype")]
    receiver.bins = rownames(receiver.df)
  } else {
    return(NULL)
  }
  ## Plotting Dimplot for the interaction:
  stRNA$Interaction_pair = "Others"
  stRNA$Interaction_pair[colnames(stRNA) %in% sender.bins] = paste0(ligand, "_", sender_Cells)
  stRNA$Interaction_pair[colnames(stRNA) %in% receiver.bins] = paste0(receptor, "_", receiver_Cells)
  plot_item = "Interaction_pair"
  prefix = paste0("BGI_interaction_CJY_dimplot/",paste0(ligand,"-",receptor),"/",id,"/",id,
                  paste0("_",sender_Cells,"-",receiver_Cells,"_distance"))
  ## check the column exist in meta.data
  if (plot_item %in% names(stRNA@meta.data)) {
    dimplot_inhouse_R(obj = stRNA, plot_item = plot_item, prefix = prefix, bg_color = bg_color)
  } else {
    tryCatch({
      print(paste0(id, " ",plot_item," not in the meta.data!"))
    }, error = function(e) {
      print(paste0("Error occurred in iteration for ", id, ": ", e$message))
    })
  }
  
  ## merge and calculate distance:
  sr.dist = rbind(sender.df, receiver.df) %>% dist() %>% dist2list() %>% .[.[,3] > 0,]
  names(sr.dist) = c("sender", "receiver", "distance")
  sr.df = as.data.frame(matrix(NA,ncol=2,nrow=0))
  names(sr.df) = c("sender_bins", "min_distance")
  for (i in unique(sr.dist$sender)){
    sam <- subset(sr.dist, sender == i)
    distance <- min(as.numeric(sam$distance))
    tmp = data.frame(sender_bins = i, min_distance = distance)
    sr.df <- rbind(sr.df, tmp)
  }
  sr.df$sender = sender_Cells
  sr.df$receiver = receiver_Cells
  return(sr.df)
}


get_surrounding_coords <- function(coord, max_distance) {
  surrounding_coords <- expand.grid(
    row = (coord$row - max_distance):(coord$row + max_distance),
    col = (coord$col - max_distance):(coord$col + max_distance)
  )
  surrounding_coords <- surrounding_coords[!((surrounding_coords$row == coord$row) &
                                               (surrounding_coords$col == coord$col)), ]
  return(surrounding_coords)
}