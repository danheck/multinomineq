
# k
#' @examples
#' debug(stratsel:::count_auto_bin)
#' cc <- stratsel:::count_auto_bin(c(k), 20, A, b, M=5000, eps = .05)
#' c(est = cc$int, sd =stratsel:::precision_count(cc$count, cc$M, log = F))
#'
#'
#' library(lineprof)
#' l <- lineprof(stratsel:::count_auto_bin(c(k), 20, A, b, eps = .05))
#' l
#' shine(l)

count_auto_bin <- function(k, n, A, b, prior = c(1, 1), M = 1000, steps = 2,
                           cmin = 10, eps = .1, maxiter = 50, start = -1){
  if (length(n) == 1)
    n <- rep(n, length(k))
  # start with many small steps and few samples
  if (length(steps) == 1)
    steps <- seq(steps, nrow(A) - 1, steps)

  # first round of samples
  if (start[1] == -1)
    start <- find_inside(A, b)
  cnt <- count_stepwise(k, n, A, b, prior, M, steps,
                        batch = M, start, progress = FALSE)
  prec <- precision_count(cnt$count, cnt$M, beta = c(.5, .5), samples = 500)

  steps <- c(steps, nrow(A))  # add last step
  iter <- 1
  while (iter < maxiter && (any(cnt$count < cmin) || prec > eps)){
    sel <- union(which(cnt$count < cmin),
                 which.min(cnt$count))
    cc <- rep(0, length(sel))

    # draw new samples (stepwise method!)
    for (i in seq_along(sel)){
      idx <- sel[i]
      if (idx > 1){
        start <- find_inside(A, b, random = TRUE)
        cc[i] <- count_step(k, n, A, b, prior, M,
                            from = steps[idx-1]-1, to = steps[idx]-1, start, FALSE)
      } else {
        ss <- seq(steps[idx])
        cc[i] <- count_binomial_cpp(k, n, A[ss,,drop = FALSE], b[ss], prior,
                                    M, batch = M + 1, progress = FALSE)["count"]
      }
    }

    # add new counts and M (independent draws!)
    cnt$count[sel] <- cnt$count[sel] + cc
    cnt$M[sel] <- cnt$M[sel] + M

    # judge accuracy
    prec <- precision_count(cnt$count, cnt$M, beta = c(.5, .5), samples = 500)
  }
  cnt
}

