
# ' Count of Samples that Satisfy Inequality Constraints
# '
# ' The class \code{ineq_count} defines a class for the Bayes factor for inequality-constrained Bayes factors.
# '
# ' @name ineq_count-class
# ' @rdname ineq_count-class
# ' @exportClass ineq_count


as_ineq_count <- function(count){

  class(count) <- c("ineq_count", "matrix")
  attr(count, "proportion") <- prod(count[,"count"]/count[,"M"])
  s_prop <- sampling_proportion(count = count[,"count"], M = count[,"M"], log = FALSE)
  attr(count, "se") <- sd(s_prop)
  count
}

#' @export
print.ineq_count <- function(x, ...){
  print.table(x)
  p <- attr(x, "proportion")
  se <- attr(x, "se")
  cat("Proporiton of samples in line with inequality constraints:\n",
      p, "( SE = ", se, ").\n")
}


