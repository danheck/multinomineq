#' Bayes Factor with Inequality and (Approximate) Equality Constraints
#'
#' To obtain the Bayes factor for the equality constraints \code{C*x = d}, a
#' sequence of approximations \code{abs(C*x - d) < delta} is used.
#'
#' @param C a matrix specifying the equality constraints \code{C*x = d} with
#'   columns refering to the free parameters (similar to \code{A})
#' @param d a vector with the same number of elements as the rows of \code{C}.
#' @param M1 number of independent samples from the encompassing model to test
#'    whether \code{A*x < b}.
#' @param M2 number of Gibbs-sampling iterations for each step of the approximation
#'   of \code{C*x = d}.
#' @param delta a vector of stepsizes that are used for the approximation.
#' @param return_Ab if \code{TRUE}, the function returns a list with the additional
#'   inequality constraints (specified via \code{A}, \code{b}, and \code{steps})
#'   that are used in the stepwise approximation \code{abs(C*x - d) < delta[i]}.
#'
#' @inheritParams sampling_multinom
#' @inheritParams bf_multinom
#'
#' @details
#' First, the encompassing Bayes factor for the equality constraint \code{A*x<b}
#' is computed using \code{M1} independent Dirichlet samples. Next, the equality
#' constraint \code{C*x=d} is approximated by drawing samples from the model
#' \code{A*x<b} and testing whether \code{abs(C*x - d) < delta[1]}. Similarly, the
#' stepsize \code{delta} is reduced step by step until \code{abs(C*x - d) < min(delta)}.
#' Klugkist  et al. (2010) show that this procedure provides a valid approximation
#' of the exact equality constraints if the step size becomes sufficiently small.
#'
#' @template ref_klugkist2010
#'
#' @examples
#' # Equality constraints:  C * x = d
#' d <- c(.5, .5, 0)
#' C <- matrix(c(1, 0, 0, 0,    # p1 = .50
#'               0, 1, 0, 0,    # p2 = .50
#'               0, 0, 1, -1),  # p3 = p4
#'             ncol = 4, byrow = TRUE)
#' k <- c(3,7, 6,4,  2,8,  5,5)
#' options <- c(2, 2, 2, 2)
#' bf_equality(k, options, C = C, d = d, delta = .5^(1:5),
#'             M1 = 50000, M2 = 5000)  # only for CRAN checks
#'
#' # check against exact equality constraints (see ?bf_binom)
#' k_binom = k[seq(1,7,2)]
#' bf_binom(k_binom, n = 10, A = matrix(0), b = 0,
#'          map = c(0, 0, 1, 1))
#' @export
bf_equality <- function(k, options, A, b, C, d, prior = rep(1, sum(options)),
                        M1 = 100000, M2 = 20000, delta = .5^(1:8),
                        return_Ab = FALSE, ...){
  if (missing(A) && missing(b)) {
    A <- matrix(NA_real_, nrow = 0, ncol = ncol(C))
    b <- numeric()
  }
  check_Abokprior(A, b, options, k, prior)
  check_Abokprior(C, d, options, k, prior)

  Ab_delta <- Ab_approx_equal(delta, C, d)
  A_delta <- rbind(A, Ab_delta$A)
  b_delta <- c(b, Ab_delta$b)
  steps <- c(nrow(A), nrow(A) + Ab_delta$steps)
  if (nrow(A) == 0)
    steps <- steps[-1]

  if (return_Ab)
    return(list(A = A_delta, b = b_delta, steps = steps))

  start <- find_inside(A_delta, b_delta)
  bf <- bf_multinom(k, options, A_delta, b_delta, prior = prior,
                    M = c(M1, M2), steps = steps, start = start, ...)
  bf["bf_00'",] <- NA  # bf_00' not defined for equality constraints
  bf
}

####################################################################
# Approximate equality constraints
#
#     |C*x - d| < delta
#
# (A)  C*x - d < delta
# (B) -C*x + d < delta
#
# (A)  C*x <  d + delta
# (B) -C*x < -d + delta

Ab_delta <- function(delta, C, d){
  A_delta <- rbind(C, -C)
  b_delta <- c(d + delta, -d + delta)
  list(A = A_delta, b = b_delta)
}

Ab_approx_equal <- function(delta, C, d){
  Ab_delta <- lapply(delta, Ab_delta, C = C, d = d)
  A <- do.call("rbind", lapply(Ab_delta, "[[", "A"))
  b <- unlist(lapply(Ab_delta, "[[", "b"))
  steps <- seq(length(delta)) * 2 * nrow(C)
  list(A = A, b = b, steps = steps)
}

####################################################################



# cnt_pr <- cnt_po <- matrix(NA, length(delta), 3,
#                            dimnames = list(NULL, c("count", "M", "steps")))
#
# cnt_pr[1,] <- count_binom(0, 0, Ab_delta[[1]]$A, Ab_delta[[1]]$b, M = M*5)
# cnt_po[1,] <- count_binom(k, n, Ab_delta[[1]]$A, Ab_delta[[1]]$b, M = M*5)
#
# for (i in seq(length(Ab_delta)-1) ){
#   prior <- sampling_binom(0, 0, Ab_delta[[i]]$A, Ab_delta[[i]]$b, M = M, progress = FALSE)
#   post  <- sampling_binom(k, n, Ab_delta[[i]]$A, Ab_delta[[i]]$b, M = M, progress = FALSE)
#
#   cnt_prior <- inside(prior, Ab_delta[[i+1]]$A, Ab_delta[[i+1]]$b)
#   cnt_post  <- inside(post,  Ab_delta[[i+1]]$A, Ab_delta[[i+1]]$b)
#
#   cnt_pr[i + 1, 1:2] <- c(sum(cnt_prior), length(cnt_prior))
#   cnt_po[i + 1, 1:2] <- c(sum(cnt_post), length(cnt_post))
#
#   print(count_to_bf(cnt_po[1:i,], cnt_pr[1:i,]))
# }
# cnt_pr[,"steps"] <- cnt_po[,"steps"] <- delta
# bf <- count_to_bf(cnt_po, cnt_pr)
