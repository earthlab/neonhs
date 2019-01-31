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

test_that('get_epsg returns the correct epsg code', {
  expect_identical(get_epsg(ex_h5), '32611')
})

test_that('get_wavelength returns the correct band names', {
  band_names <- get_wavelength(ex_h5, bands = 1:2)
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
