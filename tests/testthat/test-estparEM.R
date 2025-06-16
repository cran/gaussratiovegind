set.seed(1234)
mu1 <- 1; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x1 <- rnorm(2048, 1, 5)
y1 <- rnorm(2048, 10, 10)
z1 <- x1/y1
beta1 <- mu1/mu2; rho1 <- sigma2/sigma1; delta1 <- sigma2/mu2
theta1_0 <- estparnormratio(z1, method = "EM", eps = 1e1)
test_that("estparnormratio works 1_0", {
  expect_equal(round(theta1_0$beta), round(beta1))
  expect_equal(ceiling(theta1_0$rho), rho1)
  expect_equal(round(theta1_0$delta), delta1)
})
theta1_1 <- estparnormratio(z1, method = "EM", eps = 1e-3)
test_that("estparnormratio works 1_1", {
  expect_equal(round(theta1_1$beta, 1), beta1)
  expect_equal(round(theta1_1$rho, 0), rho1)
  expect_equal(round(theta1_1$delta, 1), delta1)
})
theta1_2 <- estparnormratio(z1, method = "EM", eps = 1e-3,
                            mux0 = 0, sigmax0 = 2, muy0 = 1, sigmay0 = 3)
test_that("estparnormratio works 1_2", {
  expect_equal(round(theta1_2$beta, 1), beta1)
  expect_equal(round(theta1_2$rho, 0), rho1)
  expect_equal(round(theta1_2$delta, 1), delta1)
})

set.seed(1234)
mu1 <- 1; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x1 <- rnorm(2048, 1, 5)
y1 <- rnorm(2048, 10, 10)
z1 <- x1/y1
beta1 <- mu1/mu2; rho1 <- sigma2/sigma1; delta1 <- sigma2/mu2
theta1_1 <- estparnormratio(z1, method = "EM", eps = 1e-3)
test_that("estparnormratio works 1_1", {
  expect_equal(round(theta1_1$beta, 1), beta1)
  expect_equal(round(theta1_1$rho, 0), rho1)
  expect_equal(round(theta1_1$delta, 1), delta1)
})
theta1_2 <- estparnormratio(z1, method = "EM", eps = 1e-3,
                            mux0 = 0, sigmax0 = 2, muy0 = 1, sigmay0 = 3)
test_that("estparnormratio works 1_2", {
  expect_equal(round(theta1_2$beta, 1), beta1)
  expect_equal(round(theta1_2$rho, 0), rho1)
  expect_equal(round(theta1_2$delta, 1), delta1)
})

set.seed(1234)
mu1 <- 0; sigma1 <- 5; mu2 <- 10; sigma2 <- 10
x2 <- rnorm(2048, mu1, sigma1)
y2 <- rnorm(2048, mu2, sigma2)
z2 <- x2/y2
beta2 <- mu1/mu2; rho2 <- sigma2/sigma1; delta2 <- sigma2/mu2
theta2_1 <- estparnormratio(z2, method = "EM", eps = 1e-3)
test_that("estparnormratio works 2_1", {
  expect_equal(round(theta2_1$beta, 1), beta2)
  expect_equal(round(theta2_1$rho, 1), rho2)
  expect_equal(round(theta2_1$delta, 1), delta2)
})
theta2_2 <- estparnormratio(z2, method = "EM", eps = 1e-3,
                            mux0 = 0, sigmax0 = 2, muy0 = 1, sigmay0 = 3)
test_that("estparnormratio works 2_2", {
  expect_equal(round(theta2_2$beta, 1), beta2)
  expect_equal(round(theta2_2$rho, 0), rho2)
  expect_equal(round(theta2_2$delta, 1), delta2)
})

set.seed(1234)
mu1 <- -1.5*20; sigma1 <- 20/3; mu2 <- 20; sigma2 <- 20/3*1.5
x3 <- rnorm(2048, mu1, sigma1)
y3 <- rnorm(2048, mu2, sigma2)
z3 <- x3/y3
beta3 <- mu1/mu2; rho3 <- sigma2/sigma1; delta3 <- sigma2/mu2
theta3_1 <- estparnormratio(z3, method = "EM", eps = 1e-4)
test_that("estparnormratio works 3_1", {
  expect_equal(round(theta3_1$beta, 1), beta3)
  expect_equal(round(theta3_1$rho, 1), rho3)
  expect_equal(round(abs(theta3_1$delta), 2), delta3)
})
theta3_2 <- estparnormratio(z3, method = "EM", eps = 1e-3,
                            mux0 = 0, sigmax0 = 2, muy0 = 1, sigmay0 = 3)
test_that("estparnormratio works 3_2", {
  expect_equal(round(theta3_2$beta, 1), beta3)
  expect_equal(round(theta3_2$rho, 1), rho3)
  expect_equal(round(abs(theta3_2$delta), 2), delta3)
})
