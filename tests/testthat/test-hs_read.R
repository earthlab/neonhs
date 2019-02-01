library(neonaop)
library(raster)
context('hs_read')

ex_h5 <- system.file('extdata', 'ex.h5', package = 'neonaop')

test_that('hs_dims returns the correct dimensions', {
  expect_identical(hs_dims(ex_h5), c(426L, 30L, 30L))
})


test_that('hs_extent returns the correct extent', {
  example_extent <- hs_extent(ex_h5)
  expected_extent <- extent(c(257000, 257030, 4111970, 4112000))
  expect_identical(example_extent, expected_extent)
})


test_that('hs_epsg returns the correct epsg code', {
  expect_identical(hs_epsg(ex_h5), '32611')
})


test_that('hs_wavelength returns the correct band names', {
  band_names <- hs_wavelength(ex_h5, bands = 1:2)
  expect_identical(band_names, c('band1_384nm', 'band2_389nm'))
})


test_that('read_hs_values provides the expected raster', {
  r <- read_hs_values(ex_h5, index = list(1:50, 1:2, 1:2))
  expect_identical(nlayers(r), 50L)
  expect_identical(xmax(r), ymax(r))
  expect_identical(xmax(r), 1)
})


test_that('calculate_index_extent generates the correct index bounds', {
  clip_ext <- extent(0, 5, 0, 5)
  h5_ext <- extent(0, 10, 0, 10)
  indices <- calculate_index_extent(extent(0, 5, 0, 5),
                                    extent(0, 10, 0, 10))
  expect_equal(indices$xmin_idx, 1)
  expect_equal(indices$xmax_idx, 5)
  expect_equal(indices$ymin_idx, 6)
  expect_equal(indices$ymax_idx, 10)
})


test_that('hs_read returns an error if given invalid bands', {
  expect_error(hs_read(ex_h5, bands = 0))
  expect_error(hs_read(ex_h5, bands = -1))
})


test_that('hs_read returns an expected raster at full extent', {
  r <- hs_read(ex_h5, bands = 1:2)
  expect_identical(raster::nlayers(r), 2L)
  expect_identical(raster::ncol(r), 30L)
  expect_identical(raster::nrow(r), 30L)
  expect_gt(min(values(r)), 0)
})


test_that('hs_read returns an expected raster at cropped extent', {
  extent_to_read <- extent(c(257000, 257020, 4111980, 4112100))
  subset_r <- hs_read(ex_h5, bands = 1:2, crop = extent_to_read)
  expect_identical(raster::nlayers(subset_r), 2L)
  expect_identical(raster::ncol(subset_r), 20L)
  expect_identical(raster::nrow(subset_r), 20L)
  expect_gt(min(values(subset_r)), 0)
})


test_that('hs_extract_pts returns NA for points outside the raster extent', {
  p <- sp::SpatialPointsDataFrame(coords = data.frame(x = 0, y = 0), 
                                  data = data.frame(id = 1),
                                  proj4string = sp::CRS('+init=epsg:32611'))
  expect_error(hs_extract_pts(ex_h5, p, bands = 1), 
               regex = 'None of the points are within')
})


test_that('hs_extract_pts returns a SpatialPointsDataFrame with right vals', {
  p <- sp::SpatialPointsDataFrame(coords = data.frame(x = c(257010, 257011),
                                                      y = c(4111990, 4111991)), 
                                  data = data.frame(id = 1:2),
                                  proj4string = sp::CRS('+init=epsg:32611'))
  band_to_read <- sample(426, 10)
  vals <- hs_extract_pts(ex_h5, p, bands = band_to_read)
  expect_identical(nrow(vals), 2L)
  expected_cols <- hs_wavelength(ex_h5, bands = band_to_read)
  expect_true(all(expected_cols %in%  names(vals)))
  expect_s4_class(vals, 'SpatialPointsDataFrame')
})
