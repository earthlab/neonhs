library(neonaop)
context('hs_read')

test_that('get_dims returns the correct dimensions', {
  ex_h5 <- system.file('extdata', 'ex.h5', package = 'neonaop')
  expect_identical(get_dims(ex_h5), c(426L, 30L, 30L))
})

