#' Data: Item Responses Theory (Karabatsos & Sheu, 2004)
#'
#' Frequency table with number of correct responses \eqn{y[ij]} on item \eqn{j}
#' obtained by the \eqn{N[ij]} respondents who have rest score \eqn{R[j] = i}.
#'
#' @format A matrix with 6 rows and columns:
#' \describe{
#'   \item{\code{rs}(i=0,...,5): }{Sum score without the item i}
#'   \item{\code{item}(j=1,...,6): }{Item number}
#' }
#' @template ref_karabatsos2004
#' @seealso The polytope for the nonparametric item response theory can be obtained using (see \code{\link{nirt_to_Ab}}).
#'
#' @examples
#' data(karabatsos2004)
#' head(karabatsos2004)
#'
"karabatsos2004"



# karabatsos2004 <- list("p" =  matrix(c(.00, .40, .40,  .00, .63,  .57,
#                                        .23, .00, .15,  .28, .58,  .57,
#                                        .23, .19, .44,  .45, .75,  .78,
#                                        .38, .47, .52,  .66, .81,  .84,
#                                        .42, .55, .52,  .89, .63,  .80,
#                                        .40, .44, .80, 1.00, .50, 1.0), 6, 6, byrow = TRUE),
#                        "N" = matrix(c(3, 5, 5, 3, 8, 7,
#                                       17, 11, 13, 18, 19, 21,
#                                       17, 21, 27, 22, 24, 23,
#                                       29, 34, 21, 35, 22, 25,
#                                       24, 20, 29, 18, 19, 20,
#                                       10, 9, 5, 4, 8, 4), 6, 6, byrow = TRUE))
# dimnames(karabatsos2004$p) <-
#   dimnames(karabatsos2004$N) <- list("rest_score" = paste0("rs", 0:5),
#                                      "item" = paste0("i", 1:6))
# karabatsos2004
# save(karabatsos2004, file = "data/karabatsos2004.RData")

# karabatsos2004$p * karabatsos2004$N
