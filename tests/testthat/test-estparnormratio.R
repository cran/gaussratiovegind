set.seed(1234)
x1 <- rnorm(2048, 1, 5)
y1 <- rnorm(2048, 10, 10)
z1 <- x1/y1
beta1 <- 0.1; rho1 <- 2; delta1 <- 1
theta1 <- estparnormratio(z1, eps = 1e-3)
test_that("estparnormratio works 1", {
  expect_equal(round(theta1$beta, 1), beta1)
  expect_equal(round(theta1$rho, 0), rho1)
  expect_equal(round(theta1$delta, 1), delta1)
})

set.seed(1234)
x2 <- rnorm(2048, 0, 5)
y2 <- rnorm(2048, 10, 10)
z2 <- x2/y2
beta2 <- 0; rho2 <- 2; delta2 <- 1
theta2 <- estparnormratio(z2, eps = 1e-3)
test_that("estparnormratio works 2", {
  expect_equal(round(theta2$beta, 1), beta2)
  expect_equal(round(theta2$rho, 1), rho2)
  expect_equal(round(theta2$delta, 1), delta2)
})

set.seed(1234)
x3 <- rnorm(2048, -1.5*20, 20/3)
y3 <- rnorm(2048, 20, 20/3*1.5)
z3 <- x3/y3
beta3 <- -1.5; rho3 <- 1.5; delta3 <- 0.5
theta3 <- estparnormratio(z3, eps = 1e-4)
test_that("estparnormratio works 3", {
  expect_equal(round(theta3$beta, 1), beta3)
  expect_equal(round(theta3$rho, 1), rho3)
  expect_equal(round(abs(theta3$delta), 2), delta3)
})
