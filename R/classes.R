#' Create a Spatial Grid Object Constructor
#'
#' @param grid_data An 'sf' spatial object containing interpolated prediction rows.
#' @param method Character string. The model engine used ("IDW" or "Regression Kriging").
#' @param target_variable Character string. The column modeled ("OC").
#' @param original_points The original 'sf' sample point dataset used to build the model.
#' @param covariate_cols Character vector. Covariates used if running Regression Kriging.
#'
#' @importFrom stats as.formula lm residuals predict
#'
#' @return A S3 object of class 'spatial_grid'.
#' @export
new_spatial_grid <- function(grid_data, method, target_variable, original_points, covariate_cols = NULL) {

  obj <- list(
    data           = grid_data,
    method         = method,
    target         = target_variable,
    points         = original_points,
    covariates     = covariate_cols,
    resolution     = 10000,
    total_cells    = nrow(grid_data)
  )

  # Assign S3 class hierarchy
  class(obj) <- c("spatial_grid", "list")
  return(obj)
}

# ===============================================

#' Print Method for Spatial Grid Objects
#'
#' @param x An object of class 'spatial_grid'.
#' @param ... Unused additional parameters.
#'
#' @export
print.spatial_grid <- function(x, ...) {
  cat("=== Geostatistical Package Grid Object ===\n")
  cat("Interpolation Method :", x$method, "\n")
  cat("Target Soil Property :", x$target, "\n")
  cat("Grid Resolution      :", x$resolution, "meters\n")
  cat("Total Active Nodes   :", x$total_cells, "cells generated\n")
  cat("------------------------------------------\n")
}

# =========================================================================

#' Summary Method for Spatial Grid Objects (Cross-Validation Evaluation)
#'
#' @param object An object of class 'spatial_grid'.
#' @param ... Unused additional parameters.
#'
#' @export
summary.spatial_grid <- function(object, ...) {

  cat("=== Running Leave-One-Out Cross-Validation Engine ===\n")
  cat("Evaluating Method:", object$method, "\n\n")

  # Routes to the specific calculation functions based on the S3 label
  if (object$method == "IDW") {

    # 'validate_idw_performance' logic
    spatial_data <- object$points
    target_col   <- object$target

    obs_coords <- sf::st_coordinates(spatial_data)
    obs_values <- spatial_data[[target_col]]
    num_points <- length(obs_values)
    cv_predictions <- numeric(num_points)

    for (i in 1:num_points) {
      target_x <- obs_coords[i, "X"]
      target_y <- obs_coords[i, "Y"]
      predictor_coords <- obs_coords[-i, , drop = FALSE]
      predictor_values <- obs_values[-i]

      dists <- sqrt((predictor_coords[, "X"] - target_x)^2 + (predictor_coords[, "Y"] - target_y)^2)
      weights <- 1 / (dists^2.0)
      cv_predictions[i] <- sum(weights * predictor_values) / sum(weights)
    }

    prediction_errors <- obs_values - cv_predictions
    mean_error <- mean(prediction_errors)
    rmse       <- sqrt(mean(prediction_errors^2))

  } else {

    # 'validate_rk_performance' logic
    spatial_data   <- object$points
    target_col     <- object$target
    covariate_cols <- object$covariates

    formula_string <- paste(target_col, "~", paste(covariate_cols, collapse = " + "))
    lm_formula     <- as.formula(formula_string)
    tabular_data   <- sf::st_drop_geometry(spatial_data)

    obs_coords <- sf::st_coordinates(spatial_data)
    obs_values <- spatial_data[[target_col]]
    num_points <- length(obs_values)
    cv_predictions <- numeric(num_points)

    for (i in 1:num_points) {
      train_tabular <- tabular_data[-i, ]
      train_coords  <- obs_coords[-i, , drop = FALSE]
      test_tabular  <- tabular_data[i, ]
      test_x        <- obs_coords[i, "X"]
      test_y        <- obs_coords[i, "Y"]

      cv_lm <- lm(lm_formula, data = train_tabular)
      train_residuals <- residuals(cv_lm)
      trend_prediction <- predict(cv_lm, newdata = test_tabular)[1]

      dists <- sqrt((train_coords[, "X"] - test_x)^2 + (train_coords[, "Y"] - test_y)^2)
      weights <- 1 / (dists^2.0)
      residual_prediction <- sum(weights * train_residuals) / sum(weights)

      cv_predictions[i] <- trend_prediction + residual_prediction
    }

    prediction_errors <- obs_values - cv_predictions
    mean_error <- mean(prediction_errors)
    rmse       <- sqrt(mean(prediction_errors^2))
  }

  # PRINT RESULTS TERMINAL SCREEN OUTPUT
  cat("--- Validation Performance Metrics ---\n")
  cat("Mean Error (ME)        :", round(mean_error, 5), "(g/kg)\n")
  cat("Root Mean Square Error :", round(rmse, 5), "(g/kg)\n")
  cat("======================================\n")
}


# =========================================================================

#' Plot Method for Spatial Grid Objects
#'
#' @param x An object of class 'spatial_grid'.
#' @param boundary_data The 'sf' spatial object used as the clipping mask (Bavaria boundary).
#' @param ... Unused additional parameters.
#'
#' @export
plot.spatial_grid <- function(x, boundary_data, ...) {

  if (!inherits(boundary_data, "sf")) {
    stop("Input 'boundary_data' must be a valid spatial 'sf' polygon object.")
  }

  # Extract the internal spatial grid data out of our S3 class wrapper
  idw_grid_output <- x$data
  polygon_mask    <- sf::st_union(boundary_data)

  if (sf::st_crs(polygon_mask) != sf::st_crs(idw_grid_output)) {
    polygon_mask <- sf::st_transform(polygon_mask, sf::st_crs(idw_grid_output))
  }

  # Clean point-in-polygon location select
  clipped_sf <- sf::st_filter(idw_grid_output, polygon_mask, .predicate = sf::st_intersects)

  # Strip geometries so geom_tile can read X and Y as clean table numbers
  flat_grid_df <- sf::st_drop_geometry(clipped_sf)

  # DYNAMICGGPLOT LAYER CONFIGURATION
  if (x$method == "IDW") {

    final_map <- ggplot2::ggplot(data = flat_grid_df, ggplot2::aes(x = X, y = Y, fill = predicted_idw)) +
      ggplot2::geom_tile(width = 10000, height = 10000) +
      ggplot2::labs(title = "Contoured Continuous IDW Interpolation Surface",
                    subtitle = "Projected Coordinates System (Meters)")
  } else {

    final_map <- ggplot2::ggplot(data = flat_grid_df, ggplot2::aes(x = X, y = Y, fill = predicted_rk_value)) +
      ggplot2::geom_tile(width = 10000, height = 10000) +
      ggplot2::labs(title = "Clipped Continuous Regression Kriging Surface",
                    subtitle = "Hybrid Environmental Trend + Spatial Residual Model")
  }

  # Add your clean cartographic map styling rules
  final_map <- final_map +
    ggplot2::scale_fill_viridis_c(option = "viridis", name = "SOC (g/kg)", na.value = NA) +
    ggplot2::geom_sf(data = polygon_mask, fill = NA, color = "black", linewidth = 0.8, inherit.aes = FALSE) +
    ggplot2::coord_sf(datum = 25832) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_text(face = "bold", size = 12),
                   panel.grid.major = ggplot2::element_line(color = "gray95"),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()
                   )

  return(final_map)
}

utils::globalVariables(c("X", "Y", "predicted_idw", "predicted_rk_value"))
#' @import dplyr
NULL
