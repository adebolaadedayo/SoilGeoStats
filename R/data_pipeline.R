
#' Load, Validate, and Clip Top Soil Data
#'
#' @param file_path Character string. Relative path to the soil CSV file.
#' @param x_col Character string. X coordinate column.
#' @param y_col Character string. Y coordinate column.
#' @param target_col Character string. The primary target variable (I have chosen "OC" in this case).
#' @param covariate_cols Character
#'  vector. Names of covariate columns.
#' @param shape_path Character string. Relative path to the boundary shapefile (.shp).
#' @return An 'sf' spatial data frame object with validated numeric columns, clipped to the shapefile if available.

load_and_validate_data <- function(
    file_path = system.file("extdata", "soil_data/LUCAS-SOIL-2018.csv", package = "SoilGeoStats"),
    x_col = "TH_LONG",
    y_col = "TH_LAT",
    target_col = "OC",
    covariate_cols = c("N", "Elev"),
    shape_path = system.file("extdata", "shape_file/Bav_Boundary.shp", package = "SoilGeoStats")
) {

  ## DATA INPUT VALIDATION: To verify argument types
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("Argument 'file_path' must be a single character string indicating a relative file path.")
  }
  if (!is.character(x_col) || !is.character(y_col) || !is.character(target_col)) {
    stop("Column name parameters must be passed strictly as character strings.")
  }

  ## FILE PATH CHECK: Ensure the CSV file exists locally
  if (!file.exists(file_path)) {
    stop(paste("File context error: No file located at relative path:", file_path))
  }

  ## DATA INGESTION: Read tabular CSV data correctly
  raw_df <- read.csv(file_path, stringsAsFactors = FALSE)

  ## STRUCTURAL VALIDATION: Check that core columns exist in the soil data
  required_cols <- c(x_col, y_col, target_col, covariate_cols)
  missing_cols <- required_cols[!required_cols %in% colnames(raw_df)]
  if (length(missing_cols) > 0) {
    stop(paste("Structural error: Missing required columns in dataset:", paste(missing_cols, collapse = ", ")))
  }

  ## DATA CLEANING: Drop any rows containing missing values (NAs) in the chosen columns

  clean_df <- raw_df

  # Loop through each column, convert to numeric.
  # If text characters exist, they will become NAs, and will be filtered subsequently.
  for (col in required_cols) {
    clean_df[[col]] <- suppressWarnings(as.numeric(raw_df[[col]]))
  }

  # Drop any rows that ended up with NAs after trying to convert text to numbers
  clean_df <- clean_df[complete.cases(clean_df[, required_cols, drop = FALSE]), ]
  if (nrow(clean_df) == 0) {
    stop("Validation failure: Zero rows remaining after removing missing data rows.")
  }

  ## GEOGRAPHIC TRANSFORMATION: Convert dataframe to an explicit sf spatial object

  # Initialize the data points using their raw WGS84 Long/Lat system (EPSG: 4326)
  spatial_sf <- sf::st_as_sf(clean_df, coords = c(x_col, y_col), crs = 4326)

  ## Transform the coordinates into the Bavarian projected system (EPSG: 25832)
  spatial_sf <- sf::st_transform(spatial_sf, crs = 25832)

  # VECTOR BOUNDARY CLIP (This is OPTIONAL - If there is a boundary shapefile)
  if (!is.null(shape_path)) {
    if (!file.exists(shape_path)) {
      stop(paste("Shapefile path error: No file located at:", shape_path))
    }

    # Read the vector boundary file
    boundary_poly <- sf::st_read(shape_path, quiet = TRUE)

    # Ensure the boundary polygon matches the coordinate reference system of the soil data points
    if (sf::st_crs(boundary_poly) != sf::st_crs(spatial_sf)) {
      boundary_poly <- sf::st_transform(boundary_poly, sf::st_crs(spatial_sf))
    }

    # Perform spatial intersection to keeps only points inside the boundary polygon
    spatial_sf <- sf::st_intersection(spatial_sf, boundary_poly)

    # Final check to ensure points actually remained after clipping
    if (nrow(spatial_sf) == 0) {
      stop("Spatial clip error: No data points fell within the provided shapefile boundary.")
    }
  }

  # RETURN VALUE: Cleaned, validated (and potentially clipped) spatial features object
  return(spatial_sf)

}
