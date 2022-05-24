#' Exponential-Family Random Graph Models for Signed Networks 
#'
#' \code{\link{sergm}}  is used to fit signed exponential-family random graph models (SERGMs), 
#' in which the probability of a given network, y, on a set of nodes is 
#' exp\{theta s(y)\}/c(theta), s(y) is a vector of network statistics for y, 
#' theta is a natural parameter vector of the same length, 
#' and c(theta) is the normalizing constant for the distribution. 
#' \code{\link{sergm}} can return a maximum pseudo-likelihood estimate or 
#' an approximate maximum likelihood estimate based on a Monte Carlo scheme. 
#' There are currently four estimation algorithms implemented, which can be chosen 
#' in the \code{\link{control.sergm}} function that has to be supplied in the \code{control} parameter. 
#' \code{MLE} is the standard MCMCLE estimation routine 
#' with a log-normal approximation as proposed by Hummel et al. for binary networks. 
#' \code{Stepping} implements the stepping algorithm also introduced by by Hummel et al., 
#' while \code{RM} evokes the Robins Monroe algorithm for signed networks. 
#' Finally, the pseudo likelihood estimates based on a multinomial approximation 
#' of the intractable likelihood can be obtain through the option \code{MPLE}. 
#' For the pseudo likelihood the standard errors are approximated by a parametric 
#' bootstrap procedure introduced by Cranmer and Desmarais for static binary ERGMs. 
#'
#' @param formula formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a numeric matrix including 
#'   the signed entries 0, 1, and -1. 
#' @param control list of parameters controling specifics of the MCMC
#' @return Object of class \code{sergm}.
#' @export
sergm = function(formula, control = control.sergm()){
  # Turn the formula into terms and network data 
  preprocess = formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network)
  # Estimation ----
  if(control$method_est == "RM"){
    # Here we initialize the coefficients 
    if(control$mple_init){
      mod_beg = mple_est(network =preprocess$network,n_actors = n_actors,mod_return = T,
                         terms = preprocess$term_names,data_list = preprocess$data_list,
                         type_list = preprocess$type_list, 
                         coef_names = preprocess$coef_names)
      coef_beg = mod_beg$coefficients
    } else {
      coef_beg = vector("numeric",length = length(preprocess$coef_names))
    }
    # Start with Phase 1 
    cat("Start with Phase 1 \n")
    phase_1 = est_var(network  = preprocess$network,
                      terms = preprocess$term_names, 
                      n_actors = n_actors,n_proposals = control$sampler_var$n_proposals,
                      n_proposals_burn_in= control$sampler_var$n_proposals_burn_in,
                      seed= control$sampler_var$seed,
                      number_networks = control$sampler_var$number_networks,
                      data_list = preprocess$data_list, type_list = preprocess$type_list,
                      coef = coef_beg,
                      mh = control$sampler_var$mh,
                      start_with_empty_net = control$sampler_var$init_empty)
    cat("Start with Phase 2 \n")
    
    D = matrix(data = 0, nrow = nrow(phase_1$var), ncol = ncol(phase_1$var))
    diag(D) = diag(phase_1$var)
    # Estimate the MLE via MCMC 
    rm_mle = rm_mle_estimation(network = preprocess$network,c = control$RM_c,
                               n_actors = n_actors,mh = control$sampler_est$mh,
                               n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                               n_proposals = control$sampler_est$n_proposals,
                               seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                               max_it = control$max_it,terms = preprocess$term_names,
                               data_list = preprocess$data_list, 
                               D = phase_1$var, 
                               tol = control$tol,type_list = preprocess$type_list,
                               beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty )
    
    
    final_coef = rm_mle[[length(rm_mle)]]
    coefficients_per_iteration = rm_mle
    
  } else if(control$method_est == "MLE"){
    # Here we initialize the coefficients 
    if(control$mple_init){
      mod_beg = mple_est(network =preprocess$network,n_actors = n_actors,mod_return = T,
                         terms = preprocess$term_names,data_list = preprocess$data_list,
                         type_list = preprocess$type_list, 
                         coef_names = preprocess$coef_names)
      coef_beg = mod_beg$coefficients
    } else {
      coef_beg = vector("numeric",length = length(preprocess$coef_names))
    }
    # Estimate the MLE via MCMC 
    cat("Start with the MLE estimation \n")
    
    mle = mle_estimation(network = preprocess$network,
                         n_actors = n_actors,mh = control$sampler_est$mh,
                         n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                         n_proposals = control$sampler_est$n_proposals,
                         seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                         max_it = control$max_it,terms = preprocess$term_names,
                         data_list = preprocess$data_list,
                         tol = control$tol,type_list = preprocess$type_list,
                         beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty )
    
    
    final_coef = mle[[length(mle)]]
    coefficients_per_iteration = mle
  } else if(control$method_est == "Stepping"){
    # Here we initialize the coefficients 
    if(control$mple_init){
      mod_beg = mple_est(network =preprocess$network,n_actors = n_actors,mod_return = T,
                         terms = preprocess$term_names,data_list = preprocess$data_list,
                         type_list = preprocess$type_list, 
                         coef_names = preprocess$coef_names)
      coef_beg = mod_beg$coefficients
    } else {
      coef_beg = vector("numeric",length = length(preprocess$coef_names))
    }
    cat("Initialize parameters with the Stepping Algorithm \n")
    # Estimate the MLE via MCMC 
    mle_stepping = mle_estimation_stepping(network = preprocess$network, steps = control$Stepping_number_grids,
                                           n_actors = n_actors,mh = control$sampler_est$mh,
                                           n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                                           n_proposals = control$sampler_est$n_proposals,
                                           seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                                           max_it = control$max_it,terms = preprocess$term_names,
                                           data_list = preprocess$data_list,
                                           tol = control$tol,type_list = preprocess$type_list,
                                           beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty )
    
    cat("Stepping Algorithm was succesfull, start MLE estimation \n")
    final_coef = mle_stepping[[length(mle_stepping)]]
    mle = mle_estimation(network = preprocess$network, 
                         n_actors = n_actors,mh = control$sampler_est$mh,
                         n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                         n_proposals = control$sampler_est$n_proposals,
                         seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                         max_it = control$max_it,terms = preprocess$term_names,
                         data_list = preprocess$data_list,
                         tol = control$tol,type_list = preprocess$type_list,
                         beg_coef = final_coef,start_with_empty_net = control$sampler_est$init_empty )
    
    final_coef = mle[[length(mle)]]
    coefficients_per_iteration = append(mle_stepping,mle)
    
  } else{
    # The last option is the MPLE
    cat("Start MPLE estimation \n")
    mod_beg = mple_est(network =preprocess$network,n_actors = n_actors,mod_return = T,
                       terms = preprocess$term_names,data_list = preprocess$data_list,
                       type_list = preprocess$type_list, 
                       coef_names = preprocess$coef_names)
    final_coef = mod_beg$coefficients
    coefficients_per_iteration = NULL
  }
  # Variance ----
  cat("Starting with the estimation of the variance\n")
  if(control$method_var == "Fisher"){
    var = est_var(network  = preprocess$network,
                  terms = preprocess$term_names, 
                  n_actors = n_actors,n_proposals = control$sampler_var$n_proposals,
                  n_proposals_burn_in= control$sampler_var$n_proposals_burn_in,
                  seed= control$sampler_var$seed,
                  number_networks = control$sampler_var$number_networks,
                  data_list = preprocess$data_list, type_list = preprocess$type_list,
                  coef = final_coef,
                  mh = control$sampler_var$mh,
                  start_with_empty_net = control$sampler_var$init_empty)
    
    names(final_coef) = preprocess$coef_names
    colnames(var$var)= preprocess$coef_names
    rownames(var$var)= preprocess$coef_names
    colnames(var$stats)= preprocess$coef_names
    mcmc_chain = mcmc(data= var$stats, start = 1, 
                      end = control$sampler_var$number_networks, thin = 1)
    t_vals = var$t_vals
    var = var$var
    
  } else if (control$method_var == "BS"){
    simulated_networks = simulate_networks_new(coefs = final_coef ,n_actors = n_actors,
                                               n_proposals_burn_in = control$sampler_var$n_proposals_burn_in,
                                               n_proposals = control$sampler_var$n_proposals,
                                               seed = control$sampler_var$seed,
                                               number_networks = control$sampler_var$number_networks,
                                               terms = preprocess$term_names,
                                               data = preprocess$data_list,
                                               type =preprocess$type_list,
                                               only_stats = FALSE,
                                               mh = control$sampler_var$mh)
    matrix_networks = lapply(simulated_networks$networks, FUN = function(x){
      map_to_mat(x$edges_pos,x$edges_neg,n_actors)
    })
    coefficients_simulated = lapply(matrix_networks, FUN = function(x){
      mple_est(network =x,n_actors = n_actors,mod_return = T,
               terms = preprocess$term_names,data_list = preprocess$data_list,
               type_list = preprocess$type_list, 
               coef_names = preprocess$coef_names)$coefficients
    })
    
    df = do.call("rbind", coefficients_simulated)
    df = sweep(df, 2, colMeans(df,na.rm = T))
    df = sweep(df, 2,final_coef,FUN = "+")
    var =  apply(df, 2, quantile, probs = c(0.025, 0.975),na.rm = T)
    coefficients_per_iteration = NULL
    t_vals = NULL
    mcmc_chain = NULL
  } else {
    if(!exists("mod_beg")){
      mod_beg = mple_est(network =preprocess$network,n_actors = n_actors,mod_return = T,
                         terms = preprocess$term_names,data_list = preprocess$data_list,
                         type_list = preprocess$type_list, 
                         coef_names = preprocess$coef_names)
    } 
    tmp_summary = summary(mod_beg)
    var = tmp_summary$cov.scaled
    coefficients_per_iteration = NULL
    t_vals = NULL
    mcmc_chain = NULL
  }
  res = list(coefficients = final_coef, 
             var = var, 
             coefficients_per_iteration = coefficients_per_iteration, 
             mcmc_chain = mcmc_chain, 
             t_vals = t_vals, 
             control = control) 
  # Evaluate the likelihood ----
  if(control$eval_likelihood){
    # browser()
    cat("Evaluate the log likelihood\n")
    sub_model = mple_est(network =preprocess$network,n_actors = n_actors,mod_return = T,
                         terms = preprocess$term_names[preprocess$dyad_idep],
                         data_list = preprocess$data_list[preprocess$dyad_idep],
                         type_list = preprocess$type_list[preprocess$dyad_idep], 
                         coef_names = preprocess$coef_names[preprocess$dyad_idep])
    final_coef_sub = final_coef
    final_coef_sub[] = 0
    # llk.dind = -sub_model$deviance/2 - -sub_model$null.deviance/2
    # llk.null = -log( (n_actors*(n_actors-1)/2)^3)
    final_coef_sub[match( names(sub_model$coefficients),names(final_coef_sub))] = sub_model$coefficients
    loglk = sergm_path_sampling(formula, coef = final_coef,coef_indep = final_coef_sub, cluster = control$cluster_bridge,
                                briges = control$n_bridges,sampler = control$sampler_path_sampling)
    # res$loglik = loglk - llk.dind + llk.null
    res$loglik = loglk + sum(log(sub_model$fitted.values[sub_model$y == 1]))
  }
  res$formula = formula
  res$n_actors = n_actors
  cat("Finished\n")
  class(res) = "sergm"
  return(res)
}

#' AIC for sergm Objects
#'
#' Prints the AIC of an estimated sergm model. 
#'
#' @param object A fitted tsergm model (object of class \code{sergm}).
#' @param ... optionally more fitted model objects (not used at the moment).
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.#' @return numeric value of AIC.
#' @export
AIC.sergm = function(object, ..., k) {
  loglik = object$loglik
  df = length(object$coefficients)
  return(k*df-2*loglik)
}



#' Model assessment for sergm Objects
#'
#' Samples from the converged object and then calculate the distribution of 
#' edgewise shared friends and enemies as well as positive and negative degree distributions. 
#'
#' @param object A fitted sergm model (object of class \code{\link{sergm}}).
#' @param sampler A fitted sergm model (object of class \code{\link{sampler.sergm}}).
#' @return numeric value of AIC.
#' @export
model_assessment = function(object, sampler = ergm.sign::sampler.sergm()){
  if(!"sergm" %in% class(object)){
    stop("Object has to be of class 'sergm'!")
  }
  preprocess = formula_preprocess(object$formula)
  # Step 1: Sample networks according to the model
  sampled_networks = sergm_simulate(object$formula,coef =object$coefficients,only_stats = F,sampler = sampler)
  # Step 2: Get the statistics of the observed network
  number_pos_edges = sum(preprocess$network == 1)
  count_ese = count_edgewise_shared_enemies(network = preprocess$network,n_actors =object$n_actors,n_edges =  number_pos_edges)
  max_ese = max(count_ese)
  count_ese = factor(count_ese, levels = 0:max_ese)
  count_esp = count_edgewise_shared_partner_pos(network = preprocess$network,n_actors =object$n_actors,n_edges =  number_pos_edges)
  max_esp = max(count_esp)
  count_esp = factor(count_esp, levels = 0:max_esp)
  degree_pos = rowSums(preprocess$network == 1)
  max_dp = max(degree_pos)
  degree_pos = factor(degree_pos, levels = 0:max_dp)
  degree_neg = rowSums(preprocess$network == -1)
  max_dn = max(degree_neg)
  degree_neg = factor(degree_neg, levels = 0:max_dn)
  # igraph_pbs = mat_to_sign_igraph(preprocess$network)
  # igraph_data_obs = c("balance_score" = balance_score(igraph_tmp),count_signed_triangles(igraph_tmp))
  
  
  # Step 3: Get statistics per sample
  statistics = lapply(sampled_networks$networks,FUN = function(x){
    network_tmp = map_to_mat(x$edges_pos,x$edges_neg,object$n_actors)
    # igraph_tmp = mat_to_sign_igraph(network_tmp)
    # igraph_data = c("balance_score" = balance_score(igraph_tmp),count_signed_triangles(igraph_tmp))
    number_pos_edges = sum(network_tmp == 1)
    count_ese = count_edgewise_shared_enemies(network = network_tmp,n_actors =object$n_actors,n_edges =  number_pos_edges)
    # max_ese = max(count_ese)
    count_esp = count_edgewise_shared_partner_pos(network = network_tmp,n_actors =object$n_actors,n_edges =  number_pos_edges)
    # max_esp = max(count_esp)
    degree_pos = rowSums(network_tmp == 1)
    # max_dp = max(degree_pos)
    degree_neg = rowSums(network_tmp == -1)
    # max_dn = max(degree_neg)
    return(list( ese = table(factor(count_ese,levels = 0:(max_ese+2))),
                 esp = table(factor(count_esp,levels = 0:(max_esp+2))),
                 dp = table(factor(degree_pos,levels = 0:(max_dp+2))),
                 dn = table(factor(degree_neg,levels = 0:(max_dn+2)))))
  })
  
  
  df_ese = do.call("rbind", lapply(statistics, FUN = function(x){data.frame(t(as.numeric(x$ese)))}))
  names(df_ese) = 1:ncol(df_ese)-1
  df_esp = do.call("rbind", lapply(statistics, FUN = function(x){data.frame(t(as.numeric(x$esp)))}))
  names(df_esp) = 1:ncol(df_esp)-1
  df_dp = do.call("rbind", lapply(statistics, FUN = function(x){data.frame(t(as.numeric(x$dp)))}))
  names(df_dp) = 1:ncol(df_dp)-1
  df_dn = do.call("rbind", lapply(statistics, FUN = function(x){data.frame(t(as.numeric(x$dn)))}))
  names(df_dn) = 1:ncol(df_dn) -1
  # df_igraph = do.call("rbind", lapply(statistics, FUN = function(x){data.frame(t(as.numeric(x$igraph_data)))}))
  # names(df_igraph) = names(statistics[[1]]$igraph_data)
  
  return(list(df_ese = df_ese,
              df_esp = df_esp,
              df_dp = df_dp,
              df_dn = df_dn,
              # df_igraph = df_igraph,
              # igraph_data_obs = igraph_data_obs,
              count_ese = table(count_ese),
              count_esp = table(count_esp),
              count_dp = table(degree_pos),
              count_dn = table(degree_neg),
              sampler = sampler))
}
#' Draw from the distribution of an Signed Exponential Family Random Graph Model
#' 
#' \code{sergm_simulate} is used to draw from exponential
#' family random network models for signed networks.  See \code{\link{sergm}} for more
#' information on these models. 
#'
#' 
#'
#' A sample of networks is randomly drawn from the specified model.  The model
#' is specified by the first argument of the function through a \code{\link{formula}}.
#' 
#' Note that the first network is sampled after \code{burnin} steps,
#' and any subsequent networks are sampled each \code{interval} steps
#' after the first.
#' 
#' 
#' @param formula A \code{\link{formula}}, which should be of the form
#' \code{y ~ <model terms>}, where \code{y} is a numeric matrix of signed entries 0, -1, and 1.  
#' For the details on the possible \code{<model terms>}, see \code{\link{sergm.terms}}.
#' @param coef Vector of parameter values for the model from which the
#'   sample is to be drawn.  
#' @param only_stats Boolean value indicating if only the matrix of the specified statistics should be returned. 
#' Otherwise, the function will also return the sampled networks represented as two lists 
#' including all negative and positives edges. This list can be transformed back to a matrix with the function \code{\link{map_to_mat}}.  
#' @param sampler An \code{\link{sampler.sergm}} object to detail how the networks are drawn, see the
#' documentation of \code{\link{sampler.sergm}} what the options are. 
#' 
#' @return Depending on the parameter \code{only_stats}, either a list including the sampled networks 
#' and a matrix of their sufficient statistics or a list only including the matrix of sufficient statistics is returned.  
#' Each network is represented as a list of two lists. One list for the negative edges and one for the positive edges. 
#' This list can be transformed back to a matrix with the function \code{\link{map_to_mat}}.   
#' 
#' @export
sergm_simulate = function(formula, coef,only_stats =T ,sampler = NULL){
  # Set default values 
  if(is.null(sampler)){
    sampler = sampler.sergm()
  }
  # Turn the formula into terms and network data 
  preprocess = formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network)
  # if(length(coef) != preprocess$term_names){
  #   stop("Number of coefficients is not equal the statistics!")
  # }
  if(sampler$init_empty == T){
    res = simulate_networks_new(coefs = coef ,n_actors = n_actors,
                                n_proposals_burn_in = sampler$n_proposals_burn_in,
                                n_proposals = sampler$n_proposals,seed = sampler$seed,
                                number_networks = sampler$number_networks,
                                terms = preprocess$term_names,
                                data = preprocess$data_list,
                                type =preprocess$type_list,only_stats = only_stats,mh = sampler$mh)
    colnames(res$stats) = preprocess$coef_names
  } else {
    res = simulation_mat(network = preprocess$network, coef = coef ,n_actors = n_actors,
                         n_proposals = sampler$n_proposals,seed = sampler$seed,
                         number_networks = sampler$number_networks,
                         terms = preprocess$term_names,
                         data_list = preprocess$data_list,
                         type_list =preprocess$type_list,
                         mh = sampler$mh)
    res = rbind(t(res[[1]]),res[[2]])
    colnames(res) = preprocess$coef_names
    res= list(stats = res)
  }
  return(res)
}

#' Summary Method for sergm Objects
#'
#' Prints basic information of the results of an estimated sergm model. 
#'
#' @param object A fitted sergm model (object of class \code{sergm}).
#' @param ... Optional parameter (currently not in use).
#' @return table of the results as a data.frame object.
#' @export
summary.sergm = function(object, ...){
  names = names(object$coefficients)
  coef = object$coefficients
  
  if(object$control$method_var == "BS"){
    ci_lower = object$var[1,]
    ci_upper = object$var[2,]
    return(data.frame(coef.names = names,
                      coef = coef,
                      ci_lower = ci_lower,
                      ci_upper = ci_upper))
  } else if(object$control$method_var == "MPLE"){
    se = sqrt(diag(object$var))
    ci_lower = coef - qnorm(p = 1-0.5*0.05)*se
    ci_upper = coef + qnorm(p = 1-0.5*0.05)*se
    return(data.frame(coef.names = names,
                      coef = coef,
                      se = se,
                      ci_lower = ci_lower,
                      ci_upper = ci_upper))
  } else{
    se = sqrt(diag(object$var))
    ci_lower = coef - qnorm(p = 1-0.5*0.05)*se
    ci_upper = coef + qnorm(p = 1-0.5*0.05)*se
    return(data.frame(coef.names = names,
                      coef = coef,
                      se = se,
                      ci_lower = ci_lower,
                      ci_upper = ci_upper))
  }
}

