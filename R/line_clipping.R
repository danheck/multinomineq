# Intersection of Line with Convex Hull
#
# p = a point within the polytope
# V = a matrix with vertices
# d = a direction vector
# dim = defines a specific direction vector along a dimension (e.g., 0,0,0,1,0,0)
#
# computes the intersection by solving for dimension-wise:
# p+lambda*d = alpha * V
line_clipping <- function (p, V, d, dim = NULL){
  if (!missing(dim) && !is.null(dim))
    d <- as.numeric(1:ncol(V) == dim)
  obj <- c(1, rep(0, nrow(V)))

  mat.pos <- rbind(rep(c(0, 1), c(1, nrow(V))), # sum(alpha) = 1
                   cbind(d, - t(V)))            # intersection p+lambda*d = sum(alpha*V)
  mat.neg <- rbind(rep(c(0, 1), c(1, nrow(V))), cbind(-d, - t(V)))
  dir <- c(rep("==", 1 + ncol(V)))
  rhs <- c(1, - p)
  lp.pos <- Rglpk_solve_LP(obj, mat.pos, dir, rhs, max = TRUE)
  lp.neg <- Rglpk_solve_LP(obj, mat.neg, dir, rhs, max = TRUE)
  # check:
  # target <- p + (lp$optimum+.000000001) * d
  # inside(target, V = V)
  rbind(p + lp.pos$optimum * d, p - lp.neg$optimum * d)
}

