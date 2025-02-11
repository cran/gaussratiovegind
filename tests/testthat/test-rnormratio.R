mux <- 1; varx <- 5
muy <- 10; vary <- 10
beta <- 0.1; rho <- 2; delta <- 1
set.seed(1234)
x <- pnorm(100000, mux, varx)
set.seed(1234)
y <- pnorm(100000, muy, vary)
z <- x/y
set.seed(1234)
zth <- rnormratio(100000, bet = beta, rho = rho, delta = delta)
test_that("dnormratio works", {
  expect_equal(ks.test(z, zth)$p.value > 0.4, TRUE)
})

mux <- 0; varx <- 5
muy <- 10; vary <- 10
beta <- 0; rho <- 2; delta <- 1
set.seed(1234)
x <- pnorm(100000, mux, varx)
set.seed(1234)
y <- pnorm(100000, muy, vary)
z <- x/y
set.seed(1234)
zth <- rnormratio(100000, bet = beta, rho = rho, delta = delta)
test_that("dnormratio2 works", {
  expect_equal(ks.test(z, zth)$p.value > 0.38, TRUE)
})
