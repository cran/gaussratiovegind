# gaussratiovegind 3.0.0

* The `dnormratio()`, `pnormratio()` and `rnormratio()` functions now work
  with the probability density function and cumulative distribution of the
  ratio of two independent or correlated Gaussian variables (in the former
  versions of gaussratiovegind, it was only independent variables).
* The `estparnormratio()` function now allows the estimation of the parameters
  of the distribution of the ratio of two independent or correlated Gaussian
  variables. A logical argument `indep` allows to indicate if the variables
  are to be considered independent or correlated.
* The `estparnormratio()`, `estparEM()` and `estparVB()` functions take new
  arguments `display`, `graph`.
* Added a `NEWS.md` file to track changes to the package.
