#' Execute Inverse Distance Weighting (IDW) Interpolation
#'
#' @param spatial_data The 'sf' spatial object containing known sample vectors.
#' @param target_col Character string. Name of the continuous variable column to model.
#' @param power Numeric. Distance decay exponent scalar weight (default for IDW = 2.0).
#' @param grid_resolution Numeric. Metric cell-grid spacing increment for mapping (meters).
#'
#' @importFrom stats as.formula lm residuals coef
#'
#' @return A data frame representing a regular grid layout with matching prediction matrix.
#' @export
interpolate_idw_surface <- function(spatial_data, target_col, power = 2.0, grid_resolution = 10000) {

  ## INPUT VALIDATION
  # Check if the input data is a valid spatial 'sf' object
  if (!inherits(spatial_data, "sf")) {
    stop("Input 'spatial_data' must be a valid spatial 'sf' object.")
  }

  # Ensure the column the user wants to map actually exists in the data layout
  if (!target_col %in% colnames(spatial_data)) {
    stop(paste("Variable mapping error: Column '", target_col, "' is missing."))
  }

  # Ensure the mathematical power parameter is a single, positive number
  if (!is.numeric(power) || power <= 0 || length(power) != 1) {
    stop("Mathematical error: IDW 'power' scalar must be a single positive numeric digit.")
  }

  # Ensure the mapping resolution is a single, positive number
  if (!is.numeric(grid_resolution) || grid_resolution <= 0 || length(grid_resolution) != 1) {
    stop("Resolution parameter must be a single positive value representing cell width.")
  }

  ## GENERATE THE EMPTY MAP GRID
  # Extract the outer boundary box coordinates (xmin, xmax, ymin, ymax) of the dataset
  bbox <- sf::st_bbox(spatial_data)

  # Create a regular sequence of coordinates from the minimum to maximum boundaries,
  # stepping forward by the grid_resolution size
  x_seq <- seq(from = bbox[["xmin"]], to = bbox[["xmax"]], by = grid_resolution)
  y_seq <- seq(from = bbox[["ymin"]], to = bbox[["ymax"]], by = grid_resolution)

  # Create an empty 2D grid matrix of XY nodes using 'expand.grid' to intersect these sequences
  prediction_grid <- expand.grid(X = x_seq, Y = y_seq)

  ## DATA PREPARATION AND STORAGE
  # Isolate the true ground coordinates (the X and Y values) of the soil data points
  obs_coords <- sf::st_coordinates(spatial_data)

  # Isolate the actual soil data values to be predicted
  obs_values <- spatial_data[[target_col]]

  # Count the total number of cells in the empty map grid
  num_cells  <- nrow(prediction_grid)

  # Create an empty numeric array to store the calculated map predictions
  predictions <- numeric(num_cells)

  ## IDW INTERPOLATION LOOP
  # Loop through every single empty coordinate cell on the grid map one by one
  for (i in 1:num_cells) {

    # Get the current grid cell's X and Y position
    gx <- prediction_grid$X[i]
    gy <- prediction_grid$Y[i]

    # Calculate the linear Euclidean distance from this current grid pixel
    # to ALL known soil sample locations at once using the Pythagorean theorem
    dists <- sqrt((obs_coords[, "X"] - gx)^2 + (obs_coords[, "Y"] - gy)^2)

    # Check if the grid pixel lands exactly on top of a real soil point.
    # If distance is 0, it can't be divided by it, so the true point value is calculated.
    if (any(dists == 0)) {
      predictions[i] <- obs_values[which(dists == 0)[1]]
      next
    }


    # Calculate the weights using the inverse distance formula: w = 1 / (d^p)
    # Points that are closer get massive weights; points that are far away get tiny weights.
    weights <- 1 / (dists^power)

    # Predicted Value = Sum of (Weights * True Values) / Sum of all Weights
    predictions[i] <- sum(weights * obs_values) / sum(weights)
  }

  ## OUTPUT
  # Append the completed vector array of calculations as a new column onto the grid
  prediction_grid$predicted_idw <- predictions

  # Make the output grid to inherit the exact spatial reference system
  grid_sf <- sf::st_as_sf(prediction_grid, coords = c("X", "Y"), crs = sf::st_crs(spatial_data), remove = FALSE)

  return(new_spatial_grid(grid_data = grid_sf, method = "IDW", target_variable = target_col, original_points = spatial_data))
}

# ==========================================================================================================================

#' Optional Multi-Variable Regression Kriging Interpolation
#'
#' @param spatial_data The 'sf' spatial object containing known soil samples.
#' @param target_col Character string. Name of dependent variable column ("OC").
#' @param covariate_cols Character vector. Covariates (e.g. c("elevation", "nitrogen")).
#' @param grid_resolution Numeric. Grid distance spacing in meters for mapping outputs.
#'
#' @importFrom stats as.formula lm residuals coef
#'
#' @return A data frame grid layout containing predicted trend components, residual surfaces, and final integrated values.
#' @export
interpolate_regression_kriging <- function(spatial_data, target_col, covariate_cols, grid_resolution = 10000) {

  ## INPUT VALIDATIloadON
  # Check that the input data is a valid spatial 'sf' object
  if (!inherits(spatial_data, "sf")) {
    stop("Input 'spatial_data' must be a valid spatial 'sf' object.")
  }

  # Verify that the primary target column ("OC") exists in the dataset
  if (!target_col %in% colnames(spatial_data)) {
    stop(paste("Target dependent variable column '", target_col, "' not found."))
  }

  # Ensures the user provided a vector of text strings identifying at least one covariate column
  if (!is.character(covariate_cols) || length(covariate_cols) < 1) {
    stop("Covariate parameters error: You must provide a character vector naming at least one predictor column.")
  }


  ## REGRESSION
  # Construct the formula string from user inputs
  # Example: "OC ~ nitrogen + elevation"
  formula_string <- paste(target_col, "~", paste(covariate_cols, collapse = " + "))

  # Convert the text string into an R statistical formula object
  lm_formula     <- as.formula(formula_string)

  # Strip away the spatial geometries temporarily to create a standard flat dataframe table.
  # This will make running a standard linear model much faster and memory-efficient.
  tabular_data <- sf::st_drop_geometry(spatial_data)

  # Fit an Ordinary Least Squares (OLS) Multiple Linear Regression mode
  regression_model <- lm(lm_formula, data = tabular_data)

  # Extract the prediction errors (residuals) for each sample point.
  # Residual = (True Observed Value) - (Regression Model Prediction).
  # These represent the local spatial variations that landscape features couldn't explain.
  spatial_data$residuals <- residuals(regression_model)

  ## GRID GENERATION
  # Find the spatial bounding coordinates of the study boundary
  bbox <- sf::st_bbox(spatial_data)

  # Generate continuous regular sequences across the mapping dimensions
  x_seq <- seq(from = bbox[["xmin"]], to = bbox[["xmax"]], by = grid_resolution)
  y_seq <- seq(from = bbox[["ymin"]], to = bbox[["ymax"]], by = grid_resolution)

  # Intersect sequences to create an empty mapping coordinate matrix grid
  prediction_grid <- expand.grid(X = x_seq, Y = y_seq)

  # Get basic array dimensions and raw matrix subsets for processing loops
  num_cells  <- nrow(prediction_grid)
  obs_coords <- sf::st_coordinates(spatial_data)
  res_values <- spatial_data$residuals

  # Create empty arrays to store intermediate and final calculations
  interpolated_residuals <- numeric(num_cells)
  fixed_trend_estimates  <- numeric(num_cells)

  # Isolate the mathematical model weights (the intercept and beta slopes coefficients)
  coefficients <- coef(regression_model)

  ## DUAL-CALCULATION LOOP (Trend + Spatial Errors)

  # Process each pixel cell on the output map grid one by one
  for (i in 1:num_cells) {
    gx <- prediction_grid$X[i]
    gy <- prediction_grid$Y[i]

    # Calculate linear Euclidean distances from this specific pixel to all known sample locations
    dists <- sqrt((obs_coords[, "X"] - gx)^2 + (obs_coords[, "Y"] - gy)^2)

    # Spatial Interpolation of Residuals (Using IDW)
    if (any(dists == 0)) {
      # If grid cell overlaps a sample perfectly, use that sample's exact error
      interpolated_residuals[i] <- res_values[which(dists == 0)[1]]
    } else {
      # Apply standard inverse distance squared weighting to smooth out the model error values
      weights <- 1 / (dists^2.0)
      interpolated_residuals[i] <- sum(weights * res_values) / sum(weights)
    }

    # Calculating the Fixed Landscape Trend Component
    # To predict the trend value at this empty grid cell, find the index of its nearest physical soil point observation. Check that point's covariate attributes (e.g., its elevation).
    closest_obs_idx <- which.min(dists)

    # Initialize the equation calculation with the static Model Intercept
    trend_val <- coefficients[1]

    # Dynamically loop through every covariate requested and apply its beta slope weight multiplication
    for (j in seq_along(covariate_cols)) {
      cov_name <- covariate_cols[j]
      cov_value <- tabular_data[closest_obs_idx, cov_name]

      # Trend Equation accumulation:
      # Trend = Intercept + (beta_1 * Covariate_1) + (beta_2 * Covariate_2)...
      trend_val <- trend_val + (coefficients[j + 1] * cov_value)
    }

    # Save the resulting landscape trend value for this pixel
    fixed_trend_estimates[i] <- trend_val
  }

  ## LAYER MERGE
  # Append the pure landscape trend surface values as a column on the grid matrix dataframe
  prediction_grid$fixed_trend <- fixed_trend_estimates

  # Append the smoothed spatial local variance error surface values as a column
  prediction_grid$interpolated_residual <- interpolated_residuals

  # Combine the environmental trend layer with the spatial error layer
  # to finalize the hybrid Regression Kriging value: Final Map = Trend + Residuals
  prediction_grid$predicted_rk_value  <- fixed_trend_estimates + interpolated_residuals

  # Return the multi-column prediction grid layout
  grid_sf <- sf::st_as_sf(prediction_grid, coords = c("X", "Y"), crs = sf::st_crs(spatial_data), remove = FALSE)

  return(new_spatial_grid(grid_data = grid_sf, method = "Regression Kriging", target_variable = target_col, original_points = spatial_data, covariate_cols = covariate_cols))
}
