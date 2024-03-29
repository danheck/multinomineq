% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_Ab_V.R
\name{Ab_multinom}
\alias{Ab_multinom}
\title{Get Constraints for Product-Multinomial Probabilities}
\usage{
Ab_multinom(options, A = NULL, b = NULL, nonneg = FALSE)
}
\arguments{
\item{options}{number of observable categories/probabilities for each item
type/multinomial distribution, e.g., \code{c(3,2)} for a ternary and binary item.}

\item{A}{a matrix defining the convex polytope via \code{A*x <= b}.
The columns of \code{A} do not include the last choice option per item type and
thus the number of columns must be equal to \code{sum(options-1)}
(e.g., the column order of \code{A} for \code{k = c(a1,a2,a2, b1,b2)}
is \code{c(a1,a2, b1)}).}

\item{b}{a vector of the same length as the number of rows of \code{A}.}

\item{nonneg}{whether to add constraints that probabilities must be nonnegative}
}
\description{
Get or add inequality constraints (or vertices) to ensure that multinomial probabilities are
positive and sum to one for all choice options within each item type.
}
\details{
If \code{A} and \code{b} are provided, the constraints are added to these inequality constraints.
}
\examples{
# three binary and two ternary choices:
options <- c(2, 2, 2, 3, 3)
Ab_multinom(options)
Ab_multinom(options, nonneg = TRUE)
}
\seealso{
\code{\link{add_fixed}}
}
