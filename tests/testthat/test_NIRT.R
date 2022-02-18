library(testthat)

mcmc.summ <- function(x)
  c(m=mean(x), sd=sd(x), quantile(x, c(.025,.975)))
################################ dichotomous

N <- 5
M <- 3

# Rasch / 1PL predictions (special case of ADISOP / ACM)
rasch <- outer(seq(-2, 2, length.out = N), -seq(-1, 2, length.out = M),
               function(x, y) (1 + exp(- x + y))^-1)

test_that("NIRT axioms match with Rasch predictions", {
  expect_silent(tmp <- nirt_to_Ab(N, M))
  A <- tmp$A
  b <- tmp$b
  expect_equal(nrow(A), length(b))
  expect_true(inside(c(rasch), A, b))

  n <- 500
  k <- rpbinom(c(rasch), n)
  # k <- rbinom(M*N, n, .5)
  start <- find_inside(A, b)
  cpost <- count_binom(k, n, A, b, M = 5000, progress = FALSE,
                       steps = seq(2,nrow(A) - 1,2), start = start)
  cprior <- count_binom(0, 0, A, b, M = 5000, progress = FALSE,
                        steps = seq(2,nrow(A) - 1,2), start = start)
  count_to_bf(cpost, cprior)
  # multinomineq:::count_auto_bin(k, n, A, b, M=5000, eps = .05)
  expect_silent(tt <- sampling_binom(k, n, A, b, M = 2000, progress =FALSE))
  plot(tt[,1], ty="l")
  expect_silent(ppp <- ppp_binom(tt, k, n))
  expect_gt(ppp[3], .01)
  #### prior distribution:
  # tt <- sampling_binom(0,0, A, b, M = 2000, prior = rep(1/(N*M), 2))
  # matrix(apply(tt, 2, mean), N)

  ################### strength of DC constraint:
  skip("functions not yet fully implemented")
  DC <- nirt_to_Ab(N, M, axioms = c("W1", "W2", "DC"))
  # count_binom(0,0, DC$exclude[[1]]$A, rep(0,3), M=5e4)
  # relative to W1 / W2
  ss <- sampling_binom(0,0, DC$A, DC$b, M = 3e4)
  dc <- sapply(DC$exclude, function(ee) inside(ss, ee$A, ee$b)) # violation DC
  table(rowSums(dc))
  mean(apply(dc, 1, any))

  spost <- sampling_binom(k, n, DC$A, DC$b, M = 3e4)
  dc.post <- sapply(DC$exclude, function(ee) inside(spost, ee$A, ee$b))
  mean(apply(dc.post, 1, any))
  spost.dc <- spost[!apply(dc.post, 1, any),]
  matrix(apply(spost.dc, 2, mean), N)
  matrix(apply(spost, 2, mean), N)
})


#############################################################
# example in Karabatsos 2001, Figure 4 p. 409:
#     Perline, Wright, and Wainer (1979), Table 2, p. 244.
ppw <- matrix(c(11,20,37,
                44,43,54,
                60,56,56), 3, 3, byrow = TRUE)
N <- matrix(c(61, 84, 82), 3, 3)
dimnames(ppw) <- dimnames(N) <- list(score=paste0("s", 3:5),
                                     item=paste0("i", c(4, 9, 2)))
p.obs <- round(ppw / N, 2)
test_that("Examples in Karabatsos (2004, Figure 4) give identical results",{
  skip("functions not yet fully implemented")
  IRT_W1 <- nirt_to_Ab(3, 3, axioms = "W1")
  pp <- sampling_binom(c(ppw), c(N), IRT_W1$A, IRT_W1$b,
                       M = 10000, burnin = 1000)
  summ <- t(apply(pp, 2, mcmc.summ))

  irt_w1 <- c(.23, .52, .73, .33, .51, .68, .55, .64, .71)
  round(cbind(p.obs = c(p.obs), irt_w1, summ), 2)
  # expect_equal(unname(summ[,1]), irt_w1,  tol = .01)   ### TODO

  ### accept-reject sampling: same results
  # M <- 2e5
  # bb <- matrix(rbeta(M*length(N), c(ppw)+1, c(N-ppw)+1), M, length(N), byrow=TRUE)
  # sel <- inside(bb, IRT_W1$A, IRT_W1$b)
  # t(apply(bb[sel,], 2, mcmc.summ))

  IRT_W2 <- nirt_to_Ab(3, 3, axioms = c("W1", "W2"))
  pp <- sampling_binom(c(ppw), c(N), IRT_W2$A, IRT_W2$b,
                       M = 10000, burnin = 1000)
  summ <- t(apply(pp, 2, mcmc.summ))
  matrix(summ[,1], 3)

  irt_w2 <- c(.18, .47, .66, .35, .55, .70, .57, .66, .74)
  round(cbind(p.obs = c(p.obs), irt_w2, summ), 2)
  # expect_equal(unname(summ[,1]), irt_w2,  tol = .01)  ### TODO
})


############# With double cancellation
# Karabatsos Figure 7, p. 415  (W1, W2, Co)
N2 <- matrix(c(84, 82, 86), 3, 3)
ppw2 <- matrix(c(15, 20, 10,
                 11, 27, 25,
                 11, 24, 55), 3, 3, byrow = TRUE)
dimnames(ppw2) <- dimnames(N2) <- list(score=paste0("s", 4:6),
                                       item=paste0("i", c(6, 1, 8)))
p.obs2 <- round(ppw2 / N2, 2)
test_that("Examples in Karabatsos (2001, Figure 7) give identical results",{
  skip("functions not yet fully implemented")
  IRT_all <- nirt_to_Ab(3, 3, axioms = c("W1", "W2", "DC"))
  pp <- sampling_binom(c(ppw2), c(N2), IRT_all$A, IRT_all$b, #prior = c(.5,.5),
                       M = 15000, burnin = 1000)

  dc.post <- sapply(IRT_all$exclude, function(ee) inside(pp, ee$A, ee$b))
  dc.violated <- apply(dc.post, 1, any)
  mean(dc.post)
  pp.dc <- pp[!dc.violated,]

  irt_all <- c(.11, .15, .21, .22, .28, .34, .33, .42, .64)
  summ <- t(apply(pp.dc, 2, mcmc.summ))
  round(cbind(p.obs2 = c(p.obs2), irt_all, summ), 2)
  # expect_equal(unname(summ[,1]), irt_all,  tol = .01)  ## TODO
})

##### karabatsos 2001, 2-PL data, p. 419, Figure 9

p <- matrix(c( .00, .28, .00, .06, .17, .50,
               .00, .36, .00, .43, .50, .71,
               .08, .46, .38, .69, .77, .62,
               .38, .50, .81, .44, .94, .94,
               .80, .45, 1.0, .80, .95, 1.0), 5, byrow = TRUE)
N <- matrix(c(18, 14,13,16,  20), 5,6)
k <- round(p*N)
test_that("Examples in Karabatsos (2001, Figure 9) give identical results",{
  irt <- nirt_to_Ab(nrow(k), ncol(k), axioms = c("W1", "W2"))
  pp <- sampling_binom(c(k), c(N), irt$A, irt$b, prior=c(1,1),
                       M = 10000, burnin = 1000)
  summ <- t(apply(pp, 2, mcmc.summ))
  summ
  matrix(summ[,1], 5)

  means <- c(.13, .64, .70, .21, .36, .85, .29, .73, .39, .80)
  subs <- c(2, 5, 10, 11, 12, 15, 16, 19, 21, 28)
  round(cbind(pp = c(p)[subs], means, summ[subs,]), 2)
  # expect_equal(unname(summ[subs,1]), means,  tol = .01)   # TODO
})


# karabatsos2004, table 3: monotonicity estimates
pmean <- matrix(c(.05, .09, .17, .07, .45, .44,
                  .19, .13, .25, .28, .57, .59,
                  .26, .23, .41, .46, .67, .73,
                  .36, .43, .50, .66, .72, .81,
                  .44, .53, .58, .86, .75, .86,
                  .55, .63, .81, .96, .80, .96), 6, byrow = TRUE)
ppp.ij <- matrix(c(.98, .09, .30, .94, .49, .60,
                   .60, .39, .56, .62, .70, .63,
                   .72, .62, .62, .62, .51, .63,
                   .66, .57, .70, .58, .42, .69,
                   .68, .69, .54, .64, .33, .50,
                   .49, .39, .74, .98, .08, .98), 6, byrow =TRUE)
test_that("Karabatsos, 2004 (Table 3) matches (posterior mean)", {
  expect_silent(data(karabatsos2004))
  IJ <- dim(karabatsos2004$k.M)
  monotonicity <- nirt_to_Ab(IJ[1], IJ[2], axioms = "W1")
  pp <- sampling_binom(k = c(karabatsos2004$k.M),
                       n = c(karabatsos2004$n.M),
                       A = monotonicity$A, b = monotonicity$b,
                       prior = c(.5, .5), M = 10000)
  post.mean <- matrix(apply(pp, 2, mean), IJ[1],
                      dimnames = dimnames(karabatsos2004$k.M))
  expect_equal(pmean, unname(post.mean), tol = .01)

  ppp <- ppp_binom(pp, c(karabatsos2004$k.M), c(karabatsos2004$n.M),
                   by = 1:prod(IJ))[,3]
  expect_equal(unname(ppp), c(ppp.ij), tol = .03)
})

# karabatsos 2004: PPP  for IIO
ppp.ij <- matrix(c(1.0, 1.0, .98, .91, .75, .45,
                   .94, .46, .72, .20, .61, .51,
                   .12, .25, .62, .62, .61, .60,
                   .59, .59, .49, .49, .59, .59,
                   .60, .36, .26, .14, .39, .58,
                   .61, .68, .44, .29, .10, .90,
                   .47, .76, .92, .98, 1.0, 1.0), 7, byrow = TRUE)
ppp.items <- c(.18, .26, .47, .21, .40, .37)
test_that("Karabatsos, 2004 (Table 6) matches (PPP values)", {
  expect_silent(data(karabatsos2004))
  IJ <- dim(karabatsos2004$k.IIO)
  monotonicity <- nirt_to_Ab(IJ[1], IJ[2], axioms = "W2")
  pp <- sampling_binom(k = c(karabatsos2004$k.IIO),
                       n = c(karabatsos2004$n.IIO),
                       A = monotonicity$A, b = monotonicity$b,
                       prior = c(.5, .5), M = 10000)
  ppp <- ppp_binom(pp, c(karabatsos2004$k.IIO), c(karabatsos2004$n.IIO),
                   by = 1:prod(IJ))
  round(matrix(ppp[,3], 7), 2)
  expect_equal(unname(ppp[,3]), c(ppp.ij), tol = .03)

  ppp <- ppp_binom(pp, c(karabatsos2004$k.IIO), c(karabatsos2004$n.IIO),
                   by = rep(1:IJ[2], each = IJ[1]))
  expect_equal(unname(ppp[,3]), c(ppp.items), tol = .03)
})

