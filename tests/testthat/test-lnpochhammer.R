test_that("pochhammer works", {
  expect_equal(lnpochhammer(0, 0), 0)
  
  expect_equal(lnpochhammer(0.5, 0), 0)
  
  expect_equal(Re(lnpochhammer(0.5, 1)), lgamma(1.5) - lgamma(0.5))
  expect_equal(Im(lnpochhammer(0.5, 1)), 0)
  
  expect_equal(Re(lnpochhammer(0.5, 2)), lgamma(2.5) - lgamma(0.5))
  expect_equal(Im(lnpochhammer(0.5, 2)), 0)
  
  expect_equal(Re(lnpochhammer(0.5, 3)), lgamma(3.5) - lgamma(0.5))
  expect_equal(Im(lnpochhammer(0.5, 3)), 0)
  
  expect_equal(lnpochhammer(7, 0), 0)
  
  expect_equal(Re(lnpochhammer(7, 1)), lgamma(8) - lgamma(7))
  expect_equal(Im(lnpochhammer(7, 1)), 0)
  
  expect_equal(Re(lnpochhammer(7, 2)), lgamma(9) - lgamma(7))
  expect_equal(Im(lnpochhammer(7, 2)), 0)
  
  expect_equal(Re(lnpochhammer(7, 3)), lgamma(10) - lgamma(7))
  expect_equal(Im(lnpochhammer(7, 3)), 0)
  
  expect_error(lnpochhammer(1, -1))
})
