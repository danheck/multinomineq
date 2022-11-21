
library(testthat)

test_that("bisection method works", {
  skip_on_cran()
  skip_if_not_installed("RcppXPtrUtils")


  # R wrapper
  bisection <- function(f, x, i, min, max, eps = 1e-10) {
    stopifnot(is.numeric(min), is.numeric(max), is.numeric(eps), i == round(i))

    if (is.function(f)) {
      stopifnot(length(f(x)) == 1, is.numeric(f(x) + 0))
      multinomineq:::bisection_r(f, x, i, min, max, eps)
    } else if (inherits(f, "XPtr")) {
      RcppXPtrUtils::checkXPtr(f, type = "SEXP", args = c("NumericVector"))
      multinomineq:::bisection_cpp(f, x, i, min, max, eps)
    } else {
      stop("'inside' must be an R function of a C++ function defined via RcppXPtrUtils::cppXPtr(...).")
    }
  }

  # define indicator function in R
  inside_r <- function(x) x[1] + x[2] - 3 < 0
  expect_equal(bisection(inside_r, c(.5, .5), 1, -10, 10), 2.5)

  # define indicator function in C++ (as pointer)
  inside_cpp <- "SEXP inside(NumericVector x){return wrap(x[0] + x[1] -3 < 0);}" # C++ indices start at 0!
  inside_ptr <- RcppXPtrUtils::cppXPtr(inside_cpp)
  expect_equal(bisection(inside_ptr, c(.5, .5), 1, -10, 10), 2.5)

  ### check visually:
  # x <- seq(-10,10,1)
  # y <- sapply(x, function(xx) inside_r(c(xx, .5)))
  # plot(x, y - .5, pch = 16) ; abline(h=0, col = 2)

  ### Comparison to simple Rcpp function: (C++ --> R --> C++)
  # Rcpp::cppFunction(inside_cpp)
  # expect_equal(bisection(inside,     c(.5,.5), 1, -10, 10), 2.5)
  # microbenchmark::microbenchmark(#cpp = bisection(inside_r,   c(.5,.5), 1, -10, 10),
  #                                r =   bisection(inside,     c(.5,.5), 1, -10, 10),
  #                                ptr = bisection(inside_ptr, c(.5,.5), 1, -10, 10))
})
