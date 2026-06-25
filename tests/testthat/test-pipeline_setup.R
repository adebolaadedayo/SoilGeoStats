test_that("Package datasets and classes initialize cleanly", {

  mock_df <- data.frame(X = c(1, 2), Y = c(3, 4), predicted_idw = c(10, 20))
  mock_sf <- sf::st_as_sf(mock_df, coords = c("X", "Y"), remove = FALSE)

  test_object <- new_spatial_grid(
    grid_data = mock_sf,
    method = "IDW",
    target_variable = "OC",
    original_points = mock_sf
  )

  expect_s3_class(test_object, "spatial_grid")
})
