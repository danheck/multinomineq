<!--
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/multinomineq)](http://cran.r-project.org/package=multinomineq)
[![Build Status](https://travis-ci.org/danheck/multinomineq.svg?branch=master)](https://travis-ci.org/danheck/multinomineq)
[![Licence](https://img.shields.io/badge/licence-GPL--2-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![monthly downloads](http://cranlogs.r-pkg.org/badges/multinomineq)](http://cranlogs.r-pkg.org/badges/multinomineq)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/multinomineq)](http://cranlogs.r-pkg.org/badges/grand-total/multinomineq)
[![Research software impact](http://depsy.org/api/package/cran/multinomineq/badge.svg)](http://depsy.org/package/r/multinomineq)
-->

R package `multinomineq`
=====

Implements Gibbs sampling and Bayes factors for multinomial models with
convex, linear-inequality constraints on the probability parameters. This
includes models that predict a linear order of binomial probabilities
(e.g., p1 < p2 < p3 < .50) and mixture models, which assume that the
parameter vector p must be inside the convex hull of a finite number of
vertices. Inequality-constrained multinomial models have applications in the
area of judgment and decision making to fit and test random utility models
(Regenwetter, M., Dana, J., & Davis-Stober, C. P. (2011).
Transitivity of preferences. Psychological Review, 118, 42–56,
<doi:10.1037/a0021150>) or to perform outcome-based strategy classification
(i.e., to select the strategy that provides the best account for a vector of
observed choice frequencies; Heck, D. W., Hilbig, B. E., & Moshagen, M.
(2017). From information processing to decisions: Formalizing and comparing
probabilistic choice models. Cognitive Psychology, 96, 26–40.
<doi:10.1016/j.cogpsych.2017.05.003>).

### Installation

To get the most recent version of `multinomineq`, the package can  directly be
installed from GitHub via:
```
# install.packages("devtools", "RcppArmadillo", "RcppProgress",
#                  "Rglpk", "quadprog")
devtools::install_github("danheck/multinomineq")
```

Note that the pacakge `rPorta` is required to transform between the vertex ($V$) and 
the inequality ($A x<b$) representation of a poyltope. The package is available on
GitHub here:  (https://github.com/TasCL/rPorta)[ https://github.com/TasCL/rPorta]

To compile C++ code, Windows and Mac require 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and 
[Xcode Command Line Tools](https://www.maketecheasier.com/install-command-line-tools-without-xcode/), respectively. 
Moreover, on Mac, it might be necessary to install the library `gfortran` manually by typing the following into the console 
([required to compile the package `RcppArmadillo`](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)):

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```
<!--The package can be downloaded from CRAN by typing `install.packages("multinomineq")` in an active R session.-->
<!--The manual is available within R by typing `vignette('multinomineq')`.-->


### Citation

If you use `multinomineq` in publications, please cite the package as follows:

Heck, D. W. (2018). multinomineq: Bayesian Inference for Inequality-Constrained
Multinomial Models. R package version 0.1.0. https://github.com/danheck/multinomineq

