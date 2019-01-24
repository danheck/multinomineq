#' @details
#' Inequality constraints are defined via an indicator function \code{inside}
#' which returns \code{inside(x)=1} (or \code{0}) if the vector of free parameters
#' \code{x} is inside (or outside) the model space. Since the vector \code{x}
#' must include only free (!) parameters, the last probability for each
#' multinomial must not be used in the function \code{inside(x)}!
#'
#' Efficiency can be improved greatly if the indicator function is defined as C++
#' code via the function \link[RcppXPtrUtils]{cppXPtr} in the package RcppXPtrUtils
#' (see below for examples). In this case, please keep in mind that indexing in C++
#' starts with 0,1,2... (not with 1,2,3,... as in R)!
