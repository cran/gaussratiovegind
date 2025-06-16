set.seed(1234)
mu1 <- 1; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x1 <- rnorm(2048, 1, 5)
y1 <- rnorm(2048, 10, 10)
z1 <- x1/y1
beta1 <- mu1/mu2; rho1 <- sigma2/sigma1; delta1 <- sigma2/mu2

theta0 <- estparnormratio(z1, method = "VB", eps = 1e1)
test_that("estparnormratio works 0", {
  expect_equal(round(theta0$beta, 1), beta1)
  expect_equal(round(theta0$rho), rho1)
  expect_equal(round(theta0$delta), delta1)
})

set.seed(1234)
mu1 <- 1; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x1 <- rnorm(2048, 1, 5)
y1 <- rnorm(2048, 10, 10)
z1 <- x1/y1
beta1 <- mu1/mu2; rho1 <- sigma2/sigma1; delta1 <- sigma2/mu2
theta1 <- estparnormratio(z1, method = "VB", eps = 1e-3)
test_that("estparnormratio works 1", {
  expect_equal(round(theta1$beta, 1), beta1)
  expect_equal(round(theta1$rho, 0), rho1)
  expect_equal(round(theta1$delta, 1), delta1)
})

set.seed(1234)
mu1 <- 0; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x2 <- rnorm(2048, mu1, sigma1)
y2 <- rnorm(2048, mu2, sigma2)
z2 <- x2/y2
beta2 <- mu1/mu2; rho2 <- sigma2/sigma1; delta2 <- sigma2/mu2
theta2 <- estparnormratio(z2, method = "VB", eps = 1e-3)
test_that("estparnormratio works 2", {
  expect_equal(round(theta2$beta, 1), beta2)
  expect_equal(round(theta2$rho, 0), rho2)
  expect_equal(round(theta2$delta, 1), delta2)
})
