
mcmc.summ <- function(x)
  c(m=mean(x), sd=sd(x), quantile(x, c(.025,.975)))
################################ dichotomous

N <- 7
M <- 5

# Rasch / 1PL predictions (special case of ADISOP / ACM)
rasch <- outer(seq(-2, 2, length.out = N), -seq(-1, 1, length.out = M),
               function(x, y) (1 + exp(- x + y))^-1)

test_that("NIRT axioms match with Rasch predictions", {
  expect_silent(tmp <- nirt_to_Ab(N, M))
  A <- tmp$A
  b <- tmp$b
  expect_equal(nrow(A), length(b))

  expect_true(inside(c(rasch), A, b))
})

#############################################################
# example in Karabatsos 2004, p. :
#     Perline, Wright, and Wainer (1979), Table 2, p. 244.
ppw <- matrix(c(11,20,37,
                44,43,54,
                60,56,56), 3, 3, byrow = TRUE)
N <- matrix(c(61, 84, 82), 3, 3)
dimnames(ppw) <- dimnames(N) <- list(score=paste0("s", 3:5),
                                     item=paste0("i", c(4, 9, 2)))
p.obs <- ppw / N

test_that("Examples in Karabatsos (2004) give identical results",{
  IRT_W1 <- nirt_to_Ab(3, 3, axioms = "WI1")
  pp <- sampling_binomial(c(ppw), c(N), IRT_W1$A, IRT_W1$b,
                          M = 10000, burnin = 1000)
  irt_w1 <- c(.23, .52, .73, .33, .51, .68, .55, .64, .71)
  summ <- t(apply(pp, 2, mcmc.summ))
  round(cbind(p.obs = c(p.obs), mean.irt_w1, summ), 2)
  expect_equal(unname(apply(pp, 2, mean)), irt_w1,  tol = .02)


  IRT_W2 <- nirt_to_Ab(3, 3, axioms = c("WI1", "WI2"))
  pp <- sampling_binomial(c(ppw), c(N), IRT_W2$A, IRT_W2$b,
                          M = 10000, burnin = 1000)
  irt_w2 <- c(.18, .47, .66, .35, .55, .70, .57, .66, .74)
  summ <- t(apply(pp, 2, mcmc.summ))
  round(cbind(p.obs = c(p.obs), irt_w2, summ), 2)
  expect_equal(unname(apply(pp, 2, mean)), irt_w2,  tol = .02)
})



k <- rpbinom(c(rasch), 20)
k <- rbinom(M*N, 20, .5)
start <- find_inside(A, b)
cpost <- count_binomial(k, 20, A, b, M = 10000, progress = FALSE,
                        steps = seq(3,99,3), start = start)
cprior <- count_binomial(0, 0, A, b, M = 10000, progress = FALSE,
                         steps = seq(3,99,3), start = start)
count_to_bf(cpost, cprior)
tt <- sampling_binomial(k, 20, A, b, M = 5000)
plot(tt[,1], ty="l")
ppp_binomial(tt, k, 20)





