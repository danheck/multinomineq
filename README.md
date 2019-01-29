<!--
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/multinomineq)](http://cran.r-project.org/package=multinomineq)
[![monthly downloads](http://cranlogs.r-pkg.org/badges/multinomineq)](http://cranlogs.r-pkg.org/badges/multinomineq)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/multinomineq)](http://cranlogs.r-pkg.org/badges/grand-total/multinomineq)
-->

[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Build Status](https://travis-ci.org/danheck/multinomineq.svg?branch=master)](https://travis-ci.org/danheck/multinomineq)
[![Coverage status](https://codecov.io/gh/danheck/multinomineq/branch/master/graph/badge.svg)](https://codecov.io/github/danheck/multinomineq?branch=master)


R package `multinomineq`
=====

<img src="man/figures/multinomineq.png" width="200" align="right">

Implements Gibbs sampling and Bayes factors for multinomial models with linear
inequality constraints on the vector of probability parameters. As special
cases, the model class includes models that predict a linear order of binomial
probabilities (e.g., p[1] < p[2] < p[3] < .50) and mixture models assuming that
the parameter vector p must be inside the convex hull of a finite number of
predicted patterns (i.e., vertices).

Inequality-constrained multinomial models have applications in the area of
judgment and decision making to fit and test random utility models (Regenwetter,
M., Dana, J., & Davis-Stober, C.P. (2011). Transitivity of preferences.
Psychological Review, 118, 42–56) or to perform outcome-based strategy
classification to select the decision strategy that provides the best account
for a vector of observed choice frequencies (Heck, D.W., Hilbig, B.E., &
Moshagen, M. (2017). From information processing to decisions: Formalizing and
comparing probabilistic choice models. Cognitive Psychology, 96, 26–40).


## References and Vignette

A formal definition of inequality-constrained multinomial models and the 
implemented computational methods for Bayesian inference is provided in:

* Heck, D. W., & Davis-Stober, C. P. (2018). 
  Multinomial models with linear inequality constraints: 
  Overview and improvements of computational methods for Bayesian inference. 
  *Manuscript under revision.* https://arxiv.org/abs/1808.07140
  
Please cite this paper if you use `multinomineq` in publications.

The package vignette provides a short introduction of how to apply the main functions of `multinomineq`:
```
vignette('multinomineq_intro')
```



## Installation

If developer tools for R are available (see below), the most recent version of 
the package `multinomineq` can directly be installed from GitHub via:
```
### install dependencies:
install.packages("devtools","RcppArmadillo","RcppProgress",
                 "Rglpk", "quadprog", "RcppXPtrUtils")

### install from Github:
devtools::install_github("danheck/multinomineq")
```

If the compilation of the source package causes any problems, the following code 
will install a binary version of `multinomineq` (only for Windows):
```
install.packages("drat")
drat::addRepo("danheck")
install.packages("multinomineq")
```

To transform between the vertex (V) and the inequality (A*x<b) representation of 
a poyltope, it is necessary to install the pacakge `rPorta`. The package is available on
GitHub (https://github.com/TasCL/rPorta) or as a precompiled package via:
```
install.packages("drat")
drat::addRepo("danheck")
install.packages("rPorta")
```


## Compilation of Source Packages

On Linux, GLPK libraries have to be installed via the console:
```
sudo apt-get install libglpk-dev
```

To compile C++ code, Windows and Mac require 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and 
[Xcode Command Line Tools](https://www.maketecheasier.com/install-command-line-tools-without-xcode/), respectively. 
Moreover, on Mac, it might be necessary to install the library `gfortran` manually by typing the following into the console 
([required to compile the package `RcppArmadillo`](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)):

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```



