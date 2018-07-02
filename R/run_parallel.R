
#' @importFrom parallel clusterExport parLapplyLB makeCluster clusterEvalQ clusterSetRNGStream
run_parallel <- function(arg, fun, cpu = 1, simplify = "count"){

  if (is.numeric(cpu)){
    ncpu <- cpu
    stop_cl <- TRUE
    cpu <- makeCluster(cpu)
    clusterSetRNGStream(cpu, sample.int(1e9, 1))
  } else {
    ncpu <- length(cpu)
    stop_cl <- FALSE
    clusterEvalQ(cpu, library("multinomineq"))
  }
  arg <- arg[!sapply(arg, is.null)]  # omit missing arguments
  arg$cpu <- 1
  clusterExport(cpu, names(arg), envir = as.environment(arg))
  out <- parLapplyLB(cpu, seq(ncpu),
                     function(i, arg) do.call(fun, arg), arg)
  if (stop_cl) stopCluster(cpu)

  if (is.null(simplify)){
    return(out)

  } else if (simplify == "count"){
    count <- do.call("++", out)
    count[,"steps"] <- out[[1]][,"steps"]
    return(count)

  } else if (simplify == "as.mcmc.list"){
    return(as.mcmc.list(out))
  } else {
    do.call(simplify, out)
  }
}


"++" <- function(...){
  if ((n <- nargs()) == 1){
    ..1
  } else {
    l <- list(...)
    do.call("++",l[-n]) + l[[n]]
  }
}
