get_surrounding_coords <- function(coord, max_distance) {
  surrounding_coords <- expand.grid(
    row = (coord$row - max_distance):(coord$row + max_distance),
    col = (coord$col - max_distance):(coord$col + max_distance)
  )
  # Apply a circular filter
  is_within_circle <- function(x, y, center_x, center_y, radius) {
    (x - center_x)^2 + (y - center_y)^2 <= radius^2
  }
  surrounding_coords <- surrounding_coords[
    is_within_circle(surrounding_coords$row, surrounding_coords$col, 
                     coord$row, coord$col, max_distance), 
  ]
  # Exclude the center coordinate
  surrounding_coords <- surrounding_coords[!((surrounding_coords$row == coord$row) &
                                               (surrounding_coords$col == coord$col)), ]
  return(surrounding_coords)
}
# Function to get surrounding coordinates for all points in the data.frame
get_all_surrounding_coords <- function(coords_df, max_distance) {
  # Initialize an empty list to store the surrounding coordinates for each point
  surrounding_coords_list <- list()
  
  # Iterate through each point in the data.frame
  for (i in 1:nrow(coords_df)) {
    coord <- list(row = coords_df[i, "row"], col = coords_df[i, "col"])
    
    # Get the surrounding coordinates for the current point
    surrounding_coords <- get_surrounding_coords(coord, max_distance)
    
    # Add the surrounding coordinates to the list
    surrounding_coords_list[[i]] <- surrounding_coords
  }
  
  # Combine all the surrounding coordinates into a single data.frame
  all_surrounding_coords <- do.call(rbind, surrounding_coords_list)
  colnames(all_surrounding_coords) <- c("row", "col")
  
  return(all_surrounding_coords)
}