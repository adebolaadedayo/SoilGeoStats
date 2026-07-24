#' Evaluate Spatial Sample Point Grid Characteristics
#'
#' @param spatial_data The 'sf' object generated from load_and_validate_data().
#' @param target_col Character string. The name of the soil variable column of interest (OC for this context).
#'
#' @return A named list containing structural summaries of geographic distances and result of exploratory data analysis.
#' @export
exploratory_data_analysis <- function(spatial_data, target_col) {

  ## DATA VALIDATION
  if (!inherits(spatial_data, "sf")) {
    stop("Method constraint: Input 'spatial_data' must be a valid spatial 'sf' object.")
  }
  if (!target_col %in% colnames(spatial_data)) {
    stop(paste("Target column '", target_col, "' does not exist within the provided data structure."))
  }

  # EXTRACT COORDINATES
  coords <- sf::st_coordinates(spatial_data)

  # DISTANCE MATRIC COMPILATION: to calculate the distance pairings across the points
  dist_matrix <- as.matrix(dist(coords))

  # EXPLORATORY DATA ANALYSIS
  soil_values   <- spatial_data[[target_col]]
  mean_metric   <- mean(soil_values) # get the mean
  median_metric <- median(soil_values) # get the median
  var_metric    <- var(soil_values) # get the variance

  # Isolate zero-point identity values to find closest pairs
  all_positive_distances <- dist_matrix[dist_matrix > 0]

  # RETURN VALUE: Data metrics arranged into a structured summary list
  summary_metrics <- list(
    total_points        = nrow(spatial_data),
    variable_mean       = mean_metric,
    variable_median     = median_metric,
    variable_variance   = var_metric,
    minimum_sampling_interval_meters = min(all_positive_distances),
    maximum_sampling_interval_meters = max(dist_matrix)
  )

  return(summary_metrics)
}
