# gaussratiovegind

## Description
It is well known that the distribution of a Gaussian ratio does not follow
a Gaussian distribution.
The lack of awareness among users of vegetation indices about this non-Gaussian
nature could lead to incorrect statistical modeling and interpretation.
This package provides tools to accurately handle and analyse such ratios:

* Probability density function
* Cumulative distribution function
* Sample simulation
* Parameter estimation

An example on the study of chlorophyll fluorescence can be found in
A. El Ghaziri et al. (2023)

## Installation

Install the package from CRAN:

```
install.packages("gaussratiovegind")
```

Or install the development version from the repository, using the [`devtools`](https://CRAN.R-project.org/package=devtools) package:

```
install.packages("devtools")
devtools::install_gitlab("imhorphen/gaussratiovegind",
                         host = "https://forgemia.inra.fr")
```

## Authors

[Pierre Santagostini](mailto:pierre.santagostini@institut-agro.fr),
[Angélina El Ghaziri](mailto:angelina.elghaziri@institut-agro.fr),
[Nizar Bouhlel](mailto:nizar.bouhlel@institut-agro.fr)
and David Rousseau

## References

El Ghaziri, A., Bouhlel, N., Sapoukhina, N., Rousseau, D.,
On the importance of non-Gaussianity in chlorophyll fluorescence imaging.
Remote Sensing 15(2), 528 (2023).
https://doi.org/10.3390/rs15020528

Bouhlel, N., Mercier, F., El Ghaziri, A., Rousseau, D., 
Parameter Estimation of the Normal Ratio Distribution with Variational Inference.
2023 31st European Signal Processing Conference (EUSIPCO), Helsinki, Finland, 2023, pp. 1823-1827.
https://doi.org/10.23919/EUSIPCO58844.2023.10290111
