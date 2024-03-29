# Gibbs Sampler for Polytope based on Vertex representation
#
# TODO: C++ implementation
#
#' @importFrom utils txtProgressBar setTxtProgressBar
sampling_V <- function(k,
                       options,
                       V,
                       prior = rep(1, sum(options)),
                       M = 5000,
                       start,
                       burnin = 10,
                       progress = TRUE) {
  stopifnot(length(M) == 1, M > 0, M == round(M), M > burnin, burnin > 0)
  options <- check_V(V, options)

  # 1) Approximate MAP estimate as starting value
  if (missing(start) || is.null(start)) {
    start <- ml_multinom(
      k = k + prior,
      V = V,
      options = options,
      n.fit = 2,
      control = list(maxit = 1000, reltol = 1e-6)
    )$p
  } else {
    stopifnot(inside(x = start, V = V))
  }

  I <- sum(options) - length(options)
  oo <- rep(1:length(options), options - 1)
  oo.all <- rep(1:length(options), options)
  shape2 <- (k + prior)[cumsum(options)]
  k_free <- drop_fixed(k, options)
  prior_free <- drop_fixed(prior, options)

  X <- matrix(NA, M, I)
  if (interactive() && progress) {
    pb <- txtProgressBar(0, M, style = 3)
  }
  X[1, ] <- start
  for (m in 2:M) {
    x <- X[m - 1, ]
    # 2) select random dimension theta[i] and get conditional distribution (dirichlet)
    idx <- sample(I)

    # 3) find lower/upper boundaries conditional on theta[-i]
    #      [see: https://en.wikipedia.org/wiki/Intersection_of_a_polyhedron_with_a_line]
    for (i in idx) {
      xi_max <- 1 - sum(x[oo == oo[i]]) + x[i]
      bnd <- line_clipping(x, V, dim = i)[, i] / xi_max
      if (abs(diff(bnd)) < 1e-15) {
        warning("Check whether parameter space/convex hull has full dimensionality (e.g., using the function: V_to_Ab)!")
      }
      x[i] <- rbeta_trunc(
        k_free[i] + prior_free[i],
        shape2[oo[i]],
        bnd[2],
        min(bnd[1], 1)
      ) * xi_max
    }
    X[m, ] <- x
    if (interactive() && progress) {
      setTxtProgressBar(pb, m)
    }
  }
  if (interactive() && progress) {
    close(pb)
  }
  colnames(X) <- colnames(V)
  X
}
