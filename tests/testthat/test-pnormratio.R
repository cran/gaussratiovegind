test_that("dnormratio works 1", {
  expect_equal(round(pnormratio(-392, 0.1, 2, 1), 3), 0)
  expect_equal(round(pnormratio(0.1, 0.1, 2, 1), 1), 0.5)
  expect_equal(round(pnormratio(390, 0.1, 2, 1), 3), 1)
})

test_that("dnormratio works 2", {
  expect_equal(round(pnormratio(-392, -1.5, 1.5, 0.5), 3), 0)
  expect_equal(round(pnormratio(-1.5, -1.5, 1.5, 0.5), 1), 0.5)
  expect_equal(round(pnormratio(384, -1.5, 1.5, 0.5), 3), 1)
})
