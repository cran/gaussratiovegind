step <- 0.001
tt <- seq(-100, 100, by = step)
par(mfrow = c(1, 2))
fth <- dnormratio(tt, 0.1, 2, 1)
test_that("dnormratio works 1", {
  expect_equal(round(sum(fth)*step, 2), 1)
})

step <- 0.001
tt <- seq(-100, 100, by = step)
par(mfrow = c(1, 2))
fth <- dnormratio(tt, -1.5, 1.5, 0.5)
pth <- cumsum(fth)*step
test_that("dnormratio works 2", {
  expect_equal(round(sum(fth)*step, 2), 1)
})
