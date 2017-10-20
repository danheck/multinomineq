
#' Gibbs Sampler for Polytope based on Vertex representation
#'
#' @export
sampling_V <- function(k, options, V, prior = rep(1, sum(options)), M = 5000,
                       start, burnin = 10, bisection = 8, progress = TRUE){

  # 1) Find inside point as starting value
  if (missing(start))
    start <- find_inside(V = V, options = options, random = TRUE)[-cumsum(options)]

  I <- sum(options) - length(options)
  X <- matrix(NA, M, I)
  oo <- rep(1:length(options), options - 1)
  oo.all <- rep(1:length(options), options)
  shape2 <- (k + prior)[cumsum(options)]
  X[1,] <- start
  for (m in 2:M){
    x <- X[m - 1,]

    # 2) select random dimension theta[i] and get conditional distribution (dirichlet)
    idx <- sample(I)

    for (i in idx){
      # 3) find lower/upper boundaries conditional on theta[-i]
      #      [see: https://en.wikipedia.org/wiki/Intersection_of_a_polyhedron_with_a_line]
      #         a) start with theta[i] = {0,1} and check with inside_V()

      ### lower/upper: search in [0, current] / [current, 1]
      l <- bisection(i, x, V, options, lower = TRUE, steps = bisection)
      u <- bisection(i, x, V, options, lower = FALSE, steps = bisection)
      ### accept-reject sampling
      x.in <- FALSE
      xi.max <- 1 - sum(x[oo == oo[i]]) + x[i]
      while(!x.in){
        # s2 <- sum((k + prior)[oo.all == oo[i]]) - (k[i] + prior[i])
        x[i] <- multinomineq:::rbeta_trunc(k[i] + prior[i], shape2[oo[i]], l, u / xi.max) * xi.max
        x.in <- inside(multinomineq:::add_fixed(x, options = options), V = V)
      }
    }
    X[m,] <- x
  }
  X
}

# find lower/upper boundary for dimension i conditional on the remaining values in x and the polytope V
bisection <- function(i, x, V, options, lower = TRUE, steps = 8){

  # maximum value on dimension i (depends on remaining probs. per option)
  oo <- rep(1:length(options), options - 1)
  xi.max <- 1 - sum(x[oo == oo[i]]) + x[i]

  x1 <- x2 <- x.test <- x
  x1.in <- x2.in <- TRUE
  if (lower){
    x1[i] <- 0
    x1.in <- inside(multinomineq:::add_fixed(x1, options = options), V = V)
  } else {
    x2[i] <-  1 * xi.max
    x2.in <- inside(multinomineq:::add_fixed(x2, options = options), V = V)
  }

  cnt <- 1
  while (cnt <= steps){
    cnt <- cnt + 1
    if (x1.in + x2.in == 1){ # exactly one is TRUE
      x.test[i] <- (x1[i] + x2[i]) / 2

      test.in <- inside(multinomineq:::add_fixed(x.test, options = options), V = V)
      if (test.in == x1.in){
        x1 <- x.test
      } else if (test.in == x2.in){
        x2 <- x.test
      } else {
        stop("logical error")
      }
      # print(x.test)
    }
  }
  ifelse(lower, min(x1[i], x2[i]), max(x1[i], x2[i]))
}

# number of "inside" evaluations:
# 2
# 2 + steps
# 2 + 2* steps

# precision of bisection:
# 1/2^(1:10)

