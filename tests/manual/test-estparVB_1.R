set.seed(1234)
mu1 <- -1.5*20; sigma1 <- 20/3; mu2 <- 20; sigma2 <- 20/3*1.5
x0 <- rnorm(2048, mu1, sigma1)
y0 <- rnorm(2048, mu2, sigma2)
z0 <- x0/y0
beta0 <- mu1/mu2; rho0 <- sigma2/sigma1; delta0 <- sigma2/mu2
theta0 <- estparnormratio(z3, method = "VB", eps = 1e-3)
test_that("estparnormratio works 0", {
  expect_equal(round(theta0$beta, 1), beta0)
  expect_equal(round(theta0$rho, 1), rho0)
  expect_equal(round(abs(theta0$delta), 1), delta0)
})

set.seed(1234)
mu1 <- 1; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x1 <- rnorm(2048, 1, 5)
y1 <- rnorm(2048, 10, 10)
z1 <- x1/y1
beta1 <- mu1/mu2; rho1 <- sigma2/sigma1; delta1 <- sigma2/mu2
theta1 <- estparnormratio(z1, method = "VB", eps = 1e-3,
                            mux0 = mu1, sigmax0 = sigma1, alphax0 = 2, betax0 = 3,
                            muy0 = mu2, sigmay0 = sigma2, alphay0 = 3, betay0 = 2)
test_that("estparnormratio works 1", {
  expect_equal(round(theta1$beta, 1), beta1)
  expect_equal(round(theta1$rho, 0), rho1)
  expect_equal(round(theta1$delta, 0), delta1)
})

set.seed(1234)
mu1 <- 0; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x2 <- rnorm(2048, mu1, sigma1)
y2 <- rnorm(2048, mu2, sigma2)
z2 <- x2/y2
beta2 <- mu1/mu2; rho2 <- sigma2/sigma1; delta2 <- sigma2/mu2
theta2 <- estparnormratio(z2, method = "VB", eps = 1e-3,
                          mux0 = mu1, sigmax0 = sigma1, alphax0 = 2, betax0 = 3,
                          muy0 = mu2, sigmay0 = sigma2, alphay0 = 3, betay0 = 2)
test_that("estparnormratio works 2", {
  expect_equal(round(theta2$beta, 1), beta2)
  expect_equal(round(theta2$rho, 0), rho2)
  expect_equal(round(theta2$delta*5, 0), delta2*5)
})

set.seed(1234)
mu1 <- -1.5*20; sigma1 <- 20/3; mu2 <- 20; sigma2 <- 20/3*1.5
x3 <- rnorm(2048, mu1, sigma1)
y3 <- rnorm(2048, mu2, sigma2)
z3 <- x3/y3
beta3 <- mu1/mu2; rho3 <- sigma2/sigma1; delta3 <- sigma2/mu2
theta3 <- estparnormratio(z3, method = "VB", eps = 1e-3,
                          mux0 = mu1, sigmax0 = sigma1, alphax0 = 2, betax0 = 3,
                          muy0 = mu2, sigmay0 = sigma2, alphay0 = 3, betay0 = 2)
test_that("estparnormratio works 3", {
  expect_equal(round(theta3$beta, 1), beta3)
  expect_equal(round(theta3$rho, 1), rho3)
  expect_equal(round(theta3$delta, 1), delta3)
})
