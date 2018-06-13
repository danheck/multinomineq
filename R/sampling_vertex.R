# Gibbs Sampler for Polytope based on Vertex representation
#
# TODO: C++ implementation
#
#' @importFrom utils txtProgressBar setTxtProgressBar
sampling_V <- function(k, options, V, prior = rep(1, sum(options)), M = 5000,
                       start, burnin = 10, progress = TRUE){
  options <- check_V(V, options)
  # 1) Approximate ML estimate as starting value
  if (missing(start) || is.null(start))
    start <- ml_multinom(k = k, V = V, options = options, n.fit = 1,
                         maxit = 50, reltol = 1e-3)$p
  else
    stopifnot(inside(x = start, V = V))

  I <- sum(options) - length(options)
  oo <- rep(1:length(options), options - 1)
  oo.all <- rep(1:length(options), options)
  shape2 <- (k + prior)[cumsum(options)]
  k_free <- drop_fixed(k, options)
  prior_free <- drop_fixed(prior, options)

  X <- matrix(NA, M, I)
  if (progress) pb <- txtProgressBar(0, M, style = 3)
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
    if (progress) setTxtProgressBar(pb, m)
  }
  if (progress) close(pb)
  colnames(X) <- colnames(V)
  X
}

