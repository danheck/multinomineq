# Gibbs Sampler for Polytope based on Vertex representation
#
# TODO: C++ implementation
#
sampling_V <- function(k, options, V, prior = rep(1, sum(options)), M = 5000,
                       start, burnin = 10, progress = TRUE){
  options <- check_V(V, options)
  # 1) Find inside point as starting value
  if (missing(start))
    start <- find_inside(V = V, options = options, random = TRUE)

  I <- sum(options) - length(options)
  oo <- rep(1:length(options), options - 1)
  oo.all <- rep(1:length(options), options)
  shape2 <- (k + prior)[cumsum(options)]
  k_free <- drop_fixed(k, options)
  prior_free <- drop_fixed(prior, options)

  X <- matrix(NA, M, I)
  X[1,] <- start
  for (m in 2:M){
    x <- X[m - 1,]
    # 2) select random dimension theta[i] and get conditional distribution (dirichlet)
    idx <- sample(I)

    # 3) find lower/upper boundaries conditional on theta[-i]
    #      [see: https://en.wikipedia.org/wiki/Intersection_of_a_polyhedron_with_a_line]
    for (i in idx){
      xi_max <- 1 - sum(x[oo == oo[i]]) + x[i]
      bnd <- line_clipping(x, V, dim = i)[,i] / xi_max
      x[i] <- rbeta_trunc(k_free[i] + prior_free[i],
                          shape2[oo[i]], bnd[2], min(bnd[1], 1)) * xi_max
    }
    X[m,] <- x
  }
  X
}
