
#' @describeIn control.sergm
#'
#' @param x See [print()].
#' @param ... Other unused parameters
#' Automatically called when an object of class \code{\link{control.sergm}} is printed.
#' @export
print.control.sergm = function(x, ...){
  object = x
  cat(paste("Estimation Method:",object$method_est," \n"))
  cat(paste("Variance Method:",object$method_var," \n"))
  cat(paste("Likelihood Evaluation:",object$eval_likelihood," \n"))
  cat(paste("MPLE Initialization:",object$mple_init," \n"))
  cat(paste("Maximal Iterations:",object$max_it," \n"))
  cat(paste("Tolerance:",object$tol," \n"))
  cat(paste("c Value in case of Robins-Monroe:",object$RM_c," \n"))
  cat(paste("Grid Number in case of Stepping Algorithm:",object$Stepping_number_grids," \n"))
  cat(paste("Number of Bridges:",object$n_bridges," \n"))
  cat(paste("Sampler for Estimation: \n"))
  print(object$sampler_est)
  
  cat(paste("\nSampler for Variance Approximation: \n"))
  print(object$sampler_var)
  cat(paste("\nSampler for Path Sampling: \n"))
  print( object$sampler_path_sampling)
  
  if(length(object$cluster) >0){
    cat(paste("\nCluster Provided: Yes\n"))
  } else {
    cat(paste("\nCluster Provided: No\n"))
  }
}
#' @describeIn sampler.sergm
#'
#' @param x See [print()].
#' @param ... Other unused parameters
#' Automatically called when an object of class \code{\link{sampler.sergm}} is printed.
#' @export
print.sampler.sergm = function(x, ...){
  object = x
  cat(paste("  Burn-In:",object$n_proposals_burn_in,"\n  MCMC Interval:",object$n_proposals,"\n  Seed:",object$seed,
            "\n  Number of Networks:",object$number_networks))
  
  cat(paste("\n  Start Empty:",object$init_empty))
  cat(paste("\n  Metropolis-Hastings:",object$mh))

}
#' @describeIn summary.tsergm
#'
#' @param x See [print()].
#' @param ... Other unused parameters
#' Automatically called when an object of class \code{\link{summary.tsergm}} is printed.
#' @export
print.summary.tsergm = function(x, ...){
  object = x
  print(object$call)
  print(object$coefficients)
}


#' @describeIn sergm
#'
#' @param x See [print()].
#' @param ... Other unused parameters
#' Automatically called when an object of class \code{\link{sergm}} is printed.
#' @export
print.sergm = function(x, ...){
  cat("SERGM\n")
  object = x
  print(object$formula)
  cat("Estimated Coefficients:")
  cat(round(as.numeric(object$coefficients),3))
  cat("\nt-Vals:",round(object$t_vals,3))
}

#' @describeIn tsergm
#'
#' @param x See [print()].
#' @param ... Other unused parameters
#' Automatically called when an object of class \code{\link{sergm}} is printed.
#' @export
print.tsergm = function(x, ...){
  cat("Temporal SERGM\n")
  object = x
  print(object$formula)
  cat("Estimated Coefficients:")
  cat(round(as.numeric(object$coefficients),3))
  cat("\nt-Vals:",round(object$t_vals,3))
}


