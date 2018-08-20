#' Data: Item Responses Theory (Karabatsos & Sheu, 2004)
#'
#' The test was part of the 1992 Trial State Assessment in Reading at
#' Grade 4, conducted by the National Assessments of Educational Progress (NAEP).
#'
#' @format A list with 4 matrices:
#' \describe{
#'   \item{\code{k.M}: }{Number of correct responses for participants with rest scores
#'       j=0,...,5 (i.e., the sum score minus the score for item i)}
#'   \item{\code{n.M}: }{Total number of participants for each cell of matrix \code{k.M}}
#'   \item{\code{k.IIO}: }{Number of correct responses for participants with sum scores j=0,...,6}
#'   \item{\code{n.IIO}: }{Total number of participants for each cell of matrix \code{k.IIO}}
#' }
#' @template ref_karabatsos2004
#' @seealso The polytope for the nonparametric item response theory can be obtained
#'     using (see \code{\link{nirt_to_Ab}}).
#'
#' @examples
#' data(karabatsos2004)
#' head(karabatsos2004)
#'
#' ######################################################
#' ##### Testing Monotonicity (M)                   #####
#' ##### (Karabatsos & Sheu, 2004, Table 3, p. 120) #####
#'
#' IJ <- dim(karabatsos2004$k.M)
#' monotonicity <- nirt_to_Ab(IJ[1], IJ[2], axioms = "W1")
#' p <- sampling_binom(k = c(karabatsos2004$k.M),
#'                     n = c(karabatsos2004$n.M),
#'                     A = monotonicity$A, b = monotonicity$b,
#'                     prior = c(.5, .5), M = 300)
#'
#' # posterior means (Table 4, p. 120)
#' post.mean <- matrix(apply(p, 2, mean), IJ[1],
#'                     dimnames = dimnames(karabatsos2004$k.M))
#' round(post.mean, 2)
#'
#' # posterior predictive checks (Table 4, p. 121)
#' ppp <- ppp_binom(p, karabatsos2004$k.M, karabatsos2004$n.M,
#'                  by = 1:prod(IJ))
#' ppp <- matrix(ppp[,3], IJ[1], dimnames = dimnames(karabatsos2004$k.M))
#' round(ppp, 2)
#'
#'
#' ######################################################
#' ##### Testing invariant item ordering (IIO)      #####
#' ##### (Karabatsos & Sheu, 2004, Table 6, p. 122) #####
#'
#' IJ <- dim(karabatsos2004$k.IIO)
#' iio <- nirt_to_Ab(IJ[1], IJ[2], axioms = "W2")
#' p <- sampling_binom(k = c(karabatsos2004$k.IIO),
#'                     n = c(karabatsos2004$n.IIO),
#'                     A = iio$A, b = iio$b,
#'                     prior = c(.5, .5), M = 300)

#' # posterior predictive checks (Table 6, p. 122)
#' ppp <- ppp_binom(prob = p, k = c(karabatsos2004$k.IIO),
#'                  n = c(karabatsos2004$n.IIO), by = 1:prod(IJ))
#' matrix(ppp[,3], 7, dimnames = dimnames(karabatsos2004$k.IIO))
#'
#' # for each item:
#' ppp <- ppp_binom(p, c(karabatsos2004$k.IIO), c(karabatsos2004$n.IIO),
#'                  by = rep(1:IJ[2], each = IJ[1]))
#' round(ppp[,3], 2)
"karabatsos2004"


######## not complete!
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
