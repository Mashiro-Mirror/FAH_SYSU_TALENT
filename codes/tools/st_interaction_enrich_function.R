gsea_term = function(gsea, filter_threshold=10, filtered_proportion = .4){
  if (nrow(gsea@result)>filter_threshold) {
    terms = gsea@result$Description[1:ceiling(nrow(gsea@result) * filtered_proportion)]
    return(terms)
  } else if (nrow(h.gsea@result) > 0) {
    terms = gsea@result$Description
    return(terms)
  } else {
    tryCatch({
      cat("has no Hallmark result!\n")
    }, error = function(e) {
      cat("Error occurred in iteration for ", id.name, ": ", e$message,"\n")
    })
  }
}

# Custom function to aggregate columns by their names
aggregate_columns_by_name <- function(df, aggregation_function = mean) {
  # Get unique column names
  unique_column_names <- unique(colnames(df))
  
  # Preallocate aggregated data frame
  aggregated_df <- matrix(0, nrow(df), length(unique_column_names))
  
  # Set column names and row names
  colnames(aggregated_df) <- unique_column_names
  rownames(aggregated_df) <- rownames(df)
  
  # Aggregate columns with the same name using the provided aggregation function
  for (column_name in unique_column_names) {
    columns_to_aggregate <- df[, colnames(df) == column_name, drop = FALSE]
    aggregated_df[, column_name] <- rowMeans(columns_to_aggregate)
  }
  
  return(aggregated_df)
}