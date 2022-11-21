# target:    bf_0u  =   marg_0 / marg_e
#   marg_0 <- .50^z * choose(n,k) * int_[Ab]  t^x * (1-t)^y dt
#   marg_e <-         choose(n,k) * int_[0,1] t^k * (1-t)^(n-k) dt
#
#
# actually computed:   bf_0u' = marg_0' / marg_e'
#   marg_0' = .50^z * choose(x+y,x) * int_[Ab]  t^x * (1-t)^y dt
#   marg_e' = .50^z * choose(x+y,x) * int_[0,1] t^x * (1-t)^y dt
#
# correction factor:
#   bf_0u = marg_0'/marg_e' *
#           .50^z * int_[0,1] t^x * (1-t)^y dt / int_[0,1] t^k * (1-t)^(n-k) dt

# map:
# 1,1,1  => same columns in A
# -2,-3  => reversed parameter (1-p) in A matrix
# 0      => constant probability .50
map_k_to_A <- function(k, n, A, map, prior = c(1, 1)) {
  if (missing(map) || is.null(map)) {
    map <- 1:ncol(A)
  }
  if (length(k) == 1) k <- rep(k, length(map))
  if (length(n) == 1) n <- rep(n, length(k))
  check_kAmap(k, A, map)

  # free parameters
  idx_free <- map == round(map) & map != 0
  free <- map
  free[!idx_free] <- NA
  k_rev <- k
  k_rev[map < 0] <- (n - k)[map < 0]
  n_aggr <- tapply(n, abs(free), sum)
  k_aggr <- tapply(k_rev, abs(free), sum)

  # not required: constants due to different binomial coefficients
  # bc_correct <- sum(lfactorial(n) - lfactorial(k) - lfactorial(n-k))
  # bc_aggr <- sum(lfactorial(n_aggr) - lfactorial(k_aggr) - lfactorial(n_aggr-k_aggr))

  # constant parameters theta_i = .50
  fix <- map == 0
  c50 <- sum(n[fix]) * log(.50)

  # constants due to different integrals over hypercubes of different dimension
  int_kn <- sum(lbeta(k + prior[1], n - k + prior[2]))
  int_aggr <- sum(lbeta(k_aggr + prior[1], n_aggr - k_aggr + prior[2]))

  list("k" = k_aggr, "n" = n_aggr, "const_map_0u" = c50 + int_aggr - int_kn)
}

get_const_map <- function(post, prior) {
  const <- 0
  if (!is.null(attr(post, "const_map_0u"))) {
    const <- attr(post, "const_map_0u")
  }
  # if (!is.null(prior$const_map_0u) && const != prior$const_map_0u)
  #   warning("Constants due to equality constraints (because of using 'map') do not match.")
  const
}

# # # # with equality constraints
# A <- matrix(c(1, -1,   # x1 < x2
#              0, 1),   # x3 < .25
#            ncol = 2, byrow = TRUE)
# b <- c(0, .25)
# # assignment to columns of A:
# map <- c(1,-1, 2,-2, .5,.5)
# k <- c(3,8, 3,6, 5,4)
# samp <- sampling_binom(k, rep(10, 6), A, b, map)
# head(samp)
# colMeans(samp)
# apply(samp, 2, plot, type = "l", ylim = c(0,.6))


check_kAmap <- function(k, A, map) {
  idx <- map[map == round(map) & map != 0] # integers
  I <- ncol(A)
  if (!all.equal(1:I, sort(unique(abs(idx))))) {
    stop("The mapping 'map' must contain each index 1,2,.., I=ncol(A) at least once.")
  }

  if (any(map != round(map))) {
    stop("Only (positive, negative, and zero) integers allowed for 'map'.")
  }
}
