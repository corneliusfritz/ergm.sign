#' Exponential-Family Random Graph Models for Dynamic Signed Networks 
#'
#' \code{\link{tsergm}}  is used to fit temporal signed exponential-family random graph models (TSERGMs), 
#' in which the probability of a given network at time point t, y_t, on a set of nodes is 
#' exp\{theta s(y_t, y_t-1)\}/c(theta, y_t-1), s(y_t, y_t-1) is a vector of network statistics for y_t, 
#' theta is a natural parameter vector of the same length, 
#' and c(theta) is the normalizing constant for the distribution. 
#' In the dynamic setting the sufficient statistics can also depend on earlier realizations of the networks 
#' and we are modelling the observed networks jointly. 
#' \code{\link{tsergm}} can return a maximum pseudo-likelihood estimate or 
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
#' @param control list of parameters controlling specifics of the MCMC (see \code{\link{control.sergm}})
#' @return Object of class \code{sergm}.
#' @export
tsergm = function(formula, control = control.sergm()){
  # Turn the formula into terms and network data 
  preprocess = t_formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network[[1]])
  # Estimation ----
  if(control$method_est == "RM"){
    # Here we initialize the coefficients 
    if(control$mple_init){
      mod_beg = t_mple_est(networks  =preprocess$network,n_actors = n_actors,mod_return = T,
                           terms = preprocess$term_names,data_list = preprocess$data_list,
                           type_list = preprocess$type_list, 
                           coef_names = preprocess$coef_names)
      coef_beg = mod_beg$coefficients
    } else {
      coef_beg = vector("numeric",length = length(preprocess$coef_names))
    }
    # Start with Phase 1 
    cat("Start with Phase 1 \n")
    
    if(length(control$cluster)>0){
      phase_1 = t_est_var_p(networks  = preprocess$network,
                            terms = preprocess$term_names, 
                            n_actors = n_actors,n_proposals = control$sampler_var$n_proposals,
                            n_proposals_burn_in= control$sampler_var$n_proposals_burn_in,
                            seed= control$sampler_var$seed,
                            number_networks = control$sampler_var$number_networks,
                            data_lists = preprocess$data_list,
                            data_lists_par = preprocess$data_list, 
                            networks_par = preprocess$network, 
                            cluster = control$cluster,
                            type_list = preprocess$type_list,
                            coef = coef_beg,
                            mh = control$sampler_var$mh,
                            start_with_empty_net = control$sampler_var$init_empty)
      
    } else{
      phase_1 = t_est_var(networks  = preprocess$network,
                          terms = preprocess$term_names, 
                          n_actors = n_actors,n_proposals = control$sampler_var$n_proposals,
                          n_proposals_burn_in= control$sampler_var$n_proposals_burn_in,
                          seed= control$sampler_var$seed,
                          number_networks = control$sampler_var$number_networks,
                          data_lists = preprocess$data_list, type_list = preprocess$type_list,
                          coef = coef_beg,
                          mh = control$sampler_var$mh,
                          start_with_empty_net = control$sampler_var$init_empty)
    } 
    
    
    cat("Start with Phase 2 \n")
    # Estimate the MLE via MCMC 
    D = matrix(data = 0, nrow = nrow(phase_1$var), ncol = ncol(phase_1$var))
    diag(D) = diag(phase_1$var)
    rm_mle = t_rm_mle_estimation(networks = preprocess$network,c = control$RM_c,
                                 n_actors = n_actors,mh = control$sampler_est$mh,
                                 n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                                 n_proposals = control$sampler_est$n_proposals,
                                 seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                                 max_it = control$max_it,terms = preprocess$term_names,
                                 data_lists = preprocess$data_list, 
                                 D =D, 
                                 tol = control$tol,type_list = preprocess$type_list,
                                 beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty )
    
    
    final_coef = rm_mle[[length(rm_mle)]]
    coefficients_per_iteration = rm_mle
    
  } else if(control$method_est == "MLE"){
    # Here we initialize the coefficients 
    if(control$mple_init){
      mod_beg = t_mple_est(networks =preprocess$network,n_actors = n_actors,mod_return = T,
                           terms = preprocess$term_names,data_lists = preprocess$data_list,
                           type_list = preprocess$type_list, 
                           coef_names = preprocess$coef_names)
      coef_beg = mod_beg$coefficients
    } else {
      coef_beg = vector("numeric",length = length(preprocess$coef_names))
    }
    # Estimate the MLE via MCMC 
    cat("Start with the MLE estimation \n")
    
    if(length(control$cluster)>0){
      mle = t_mle_estimation_p(networks = preprocess$network,
                               n_actors = n_actors,mh = control$sampler_est$mh,
                               n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                               n_proposals = control$sampler_est$n_proposals,
                               seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                               max_it = control$max_it,terms = preprocess$term_names,
                               data_lists = preprocess$data_list,
                               cluster = control$cluster,
                               data_lists_par = preprocess$data_list,
                               networks_par = preprocess$network,
                               tol = control$tol,type_list = preprocess$type_list,
                               beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty )
    } 
    else {
      mle = t_mle_estimation(networks = preprocess$network,
                             n_actors = n_actors,mh = control$sampler_est$mh,
                             n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                             n_proposals = control$sampler_est$n_proposals,
                             seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                             max_it = control$max_it,terms = preprocess$term_names,
                             data_lists = preprocess$data_list,
                             tol = control$tol,type_list = preprocess$type_list,
                             beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty )
    }
    
    final_coef = mle[[length(mle)]]
    coefficients_per_iteration = mle
  } else if(control$method_est == "Stepping"){
    # Here we initialize the coefficients 
    if(control$mple_init){
      mod_beg = t_mple_est(networks =preprocess$network,n_actors = n_actors,mod_return = T,
                           terms = preprocess$term_names,data_lists = preprocess$data_list,
                           type_list = preprocess$type_list, 
                           coef_names = preprocess$coef_names)
      coef_beg = mod_beg$coefficients
    } else {
      coef_beg = vector("numeric",length = length(preprocess$coef_names))
    }
    cat("Initialize parameters with the Stepping Algorithm \n")
    # Estimate the MLE via MCMC 
    if(length(control$cluster)>0){
      mle_stepping = t_mle_estimation_stepping_p(networks = preprocess$network,
                                                 steps = control$Stepping_number_grids,
                                                 n_actors = n_actors,
                                                 mh = control$sampler_est$mh,
                                                 n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                                                 n_proposals = control$sampler_est$n_proposals,
                                                 seed = control$sampler_est$seed,
                                                 number_networks = control$sampler_est$number_networks,
                                                 max_it = control$max_it,terms = preprocess$term_names,
                                                 data_lists = preprocess$data_list,
                                                 data_lists_par = preprocess$data_list,
                                                 networks_par = preprocess$network,
                                                 tol = control$tol,
                                                 type_list = preprocess$type_list,
                                                 beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty,
                                                 cluster =  control$cluster)
    } else {
      mle_stepping = t_mle_estimation_stepping(networks = preprocess$network, steps = control$Stepping_number_grids,
                                               n_actors = n_actors,mh = control$sampler_est$mh,
                                               n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                                               n_proposals = control$sampler_est$n_proposals,
                                               seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                                               max_it = control$max_it,terms = preprocess$term_names,
                                               data_lists = preprocess$data_list,
                                               tol = control$tol,type_list = preprocess$type_list,
                                               beg_coef = coef_beg,start_with_empty_net = control$sampler_est$init_empty)
      
    }
    
    cat("Stepping Algorithm was succesfull, start MLE estimation \n")
    final_coef = mle_stepping[[length(mle_stepping)]]
    if(length(control$cluster)>0){
      mle = t_mle_estimation_p(networks = preprocess$network,
                               n_actors = n_actors,mh = control$sampler_est$mh,
                               n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                               n_proposals = control$sampler_est$n_proposals,
                               seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                               max_it = control$max_it,terms = preprocess$term_names,
                               data_lists = preprocess$data_list,
                               networks_par = preprocess$network, 
                               cluster = control$cluster,
                               data_lists_par = preprocess$data_list,
                               tol = control$tol,type_list = preprocess$type_list,
                               beg_coef = final_coef,start_with_empty_net = control$sampler_est$init_empty )
    } 
    else {
      mle = t_mle_estimation(networks = preprocess$network,
                             n_actors = n_actors,mh = control$sampler_est$mh,
                             n_proposals_burn_in = control$sampler_est$n_proposals_burn_in,
                             n_proposals = control$sampler_est$n_proposals,
                             seed = control$sampler_est$seed,number_networks = control$sampler_est$number_networks,
                             max_it = control$max_it,terms = preprocess$term_names,
                             data_lists = preprocess$data_list,
                             tol = control$tol,type_list = preprocess$type_list,
                             beg_coef = final_coef,start_with_empty_net = control$sampler_est$init_empty )
    }
    final_coef = mle[[length(mle)]]
    coefficients_per_iteration = append(mle_stepping,mle)
    
  } else{
    # The last option is the MPLE
    cat("Start MPLE estimation \n")
    mod_beg = t_mple_est(networks = preprocess$network,n_actors = n_actors,mod_return = T,
                         terms = preprocess$term_names,data_lists = preprocess$data_list,
                         type_list = preprocess$type_list, 
                         coef_names = preprocess$coef_names)
    final_coef = mod_beg$coefficients
    coefficients_per_iteration = NULL
  }
  # Variance ----
  cat("Starting with the estimation of the variance\n")
  if(control$method_var == "Fisher"){
    if(length(control$cluster)>0){
      var = t_est_var_p(networks  = preprocess$network,
                        terms = preprocess$term_names, 
                        n_actors = n_actors,n_proposals = control$sampler_var$n_proposals,
                        n_proposals_burn_in= control$sampler_var$n_proposals_burn_in,
                        seed= control$sampler_var$seed,
                        cluster = control$cluster,
                        data_lists_par = preprocess$data_list,
                        number_networks = control$sampler_var$number_networks,
                        data_lists = preprocess$data_list, 
                        networks_par = preprocess$network, 
                        coef = final_coef,
                        type_list = preprocess$type_list,
                        mh = control$sampler_var$mh,
                        start_with_empty_net = control$sampler_var$init_empty)
    } else {
      var = t_est_var(networks  = preprocess$network,
                      terms = preprocess$term_names, 
                      n_actors = n_actors,n_proposals = control$sampler_var$n_proposals,
                      n_proposals_burn_in= control$sampler_var$n_proposals_burn_in,
                      seed= control$sampler_var$seed,
                      number_networks = control$sampler_var$number_networks,
                      data_lists = preprocess$data_list, type_list = preprocess$type_list,
                      coef = final_coef,
                      mh = control$sampler_var$mh,
                      start_with_empty_net = control$sampler_var$init_empty)
    }
    names(final_coef) = preprocess$coef_names
    colnames(var$var)= preprocess$coef_names
    rownames(var$var)= preprocess$coef_names
    colnames(var$stats_mean)= preprocess$coef_names
    colnames(var$stats_raw)= preprocess$coef_names
    
    mcmc_chain = mcmc(data= var$stats_mean, start = 1, 
                      end = control$sampler_var$number_networks, thin = 1)
    mcmc_chain_raw = mcmc(data= var$stats_raw, start = 1, 
                          end = control$sampler_var$number_networks, thin = 1)
    t_vals = var$t_vals
    var = var$var
    # browser()
    
    # Calculate the MCMC error
    obs_stat = colSums(temporal_count_statistics(formula))
    coefficients_0 = final_coef
    coefficients  = final_coef + solve(var(mcmc_chain),obs_stat- colMeans(mcmc_chain_raw))
    coefficients_per_iteration = append(coefficients_per_iteration,final_coef)
    
    final_coef = coefficients
    mcmc_var = sergm.MCMCse(theta0 = coefficients_0,
                            theta = coefficients,
                            statsmatrices = mcmc_chain,
                            H = var(mcmc_chain))
    
    # M = nrow(mcmc_chain_raw)
    # autocorr_vals = autocorr((obs_stat-mcmc_chain_raw)*exp((coefficients - coefficients_0)*sweep(mcmc_chain, 2, apply(mcmc_chain,2,max))),lags = 0:30)
    # V = 1/(M^2)*colSums(exp((coefficients_0- coefficients)*mcmc_chain_raw))*apply(autocorr_vals,c(2,3),sum)*2
    # mcmc_var = 1/M*var %*%V %*%var
    colnames(mcmc_var)= preprocess$coef_names
    rownames(mcmc_var)= preprocess$coef_names
    # 
    
  } else if (control$method_var == "BS"){
    # TODO: This needs to be done for each network 
    if(length(control$cluster)>0){
      simulated_networks = parLapply(cl =control$cluster,X = 1:length(preprocess$data_list),
                                     fun = function(x){
                                       simulate_networks_new(coefs = final_coef ,n_actors = n_actors,
                                                             n_proposals_burn_in = control$sampler_var$n_proposals_burn_in,
                                                             n_proposals = control$sampler_var$n_proposals,
                                                             seed = control$sampler_var$seed,
                                                             number_networks = control$sampler_var$number_networks,
                                                             terms = preprocess$term_names,
                                                             data = preprocess$data_list[[x]],
                                                             type =preprocess$type_list,
                                                             only_stats = FALSE,
                                                             mh = control$sampler_var$mh)
                                     })
    } else{
      simulated_networks = lapply(X = 1:length(preprocess$data_list),
                                  FUN = function(x){
                                    simulate_networks_new(coefs = final_coef ,n_actors = n_actors,
                                                          n_proposals_burn_in = control$sampler_var$n_proposals_burn_in,
                                                          n_proposals = control$sampler_var$n_proposals,
                                                          seed = control$sampler_var$seed,
                                                          number_networks = control$sampler_var$number_networks,
                                                          terms = preprocess$term_names,
                                                          data = preprocess$data_list[[x]],
                                                          type =preprocess$type_list,
                                                          only_stats = FALSE,
                                                          mh = control$sampler_var$mh)
                                  })
    }
    matrix_networks = lapply(simulated_networks,FUN = function(y){
      lapply(y$networks, FUN = function(x){
        map_to_mat(x$edges_pos,x$edges_neg,n_actors)
      })
    })
    matrix_networks = do.call(rbind, matrix_networks)
    if(length(control$cluster)>0){
      coefficients_simulated = parApply(cl = control$cluster,
                                        matrix_networks,MARGIN = 2, FUN = function(x){
                                          t_mple_est(networks =x,n_actors = n_actors,mod_return = T,
                                                     terms = preprocess$term_names,data_lists = preprocess$data_list,
                                                     type_list = preprocess$type_list, 
                                                     coef_names = preprocess$coef_names)$coefficients
                                        })
    } else {
      coefficients_simulated = apply(matrix_networks,MARGIN = 2, FUN = function(x){
        t_mple_est(networks =x,n_actors = n_actors,mod_return = T,
                   terms = preprocess$term_names,data_lists = preprocess$data_list,
                   type_list = preprocess$type_list, 
                   coef_names = preprocess$coef_names)$coefficients
      })
    }
    
    # browser()
    df = data.frame(t(coefficients_simulated))
    # Subtract the mean -> center
    df = sweep(df, 2, colMeans(df,na.rm = T))
    # Average around the final coefficient
    df = sweep(df, 2,final_coef,FUN = "+")
    var =  apply(df, 2, quantile, probs = c(0.025, 0.975),na.rm = T)
    coefficients_per_iteration = NULL
    t_vals = NULL
    mcmc_chain = NULL
    mcmc_chain_raw = NULL
    mcmc_var = NULL
  } else {
    if(!exists("mod_beg")){
      mod_beg = t_mple_est(networks =preprocess$network,n_actors = n_actors,mod_return = T,
                           terms = preprocess$term_names,data_lists = preprocess$data_list,
                           type_list = preprocess$type_list, 
                           coef_names = preprocess$coef_names)
    } 
    tmp_summary = summary(mod_beg)
    var = tmp_summary$cov.scaled
    coefficients_per_iteration = NULL
    t_vals = NULL
    mcmc_chain = NULL
    mcmc_chain_raw = NULL
    mcmc_var = NULL
    
  }
  res = list(coefficients = final_coef, 
             var = var, 
             mcmc_var = mcmc_var,
             coefficients_per_iteration = coefficients_per_iteration, 
             mcmc_chain = mcmc_chain, 
             mcmc_chain_raw = mcmc_chain_raw,
             t_vals = t_vals, 
             control = control, 
             call = sys.calls()[[1]]) 
  # Evaluate the likelihood ----
  if(control$eval_likelihood){
    # browser()
    cat("Evaluate the log likelihood\n")
    data_list_dyad = list()
    for(i in 1:length(preprocess$network)){
      data_list_dyad[[i]] = preprocess$data_list[[i]][preprocess$dyad_idep]
    }
    
    sub_model = t_mple_est(networks =preprocess$network,n_actors = n_actors,mod_return = T,
                           terms = preprocess$term_names[preprocess$dyad_idep],
                           data_lists = data_list_dyad,
                           type_list = preprocess$type_list[preprocess$dyad_idep], 
                           coef_names = preprocess$coef_names[preprocess$dyad_idep])
    final_coef_sub = final_coef
    final_coef_sub[] = 0
    # llk.dind = -sub_model$deviance/2 - -sub_model$null.deviance/2
    # llk.null = -log( (n_actors*(n_actors-1)/2)^3)
    final_coef_sub[match( names(sub_model$coefficients),names(final_coef_sub))] = sub_model$coefficients
    loglk_ls = list()
    terms = attr(terms.formula(formula),"term.labels")
    indicator_data = grep(pattern = "data =",x = terms)
    formulae = preprocess_list = global_stats = list()
    for(i in 1:length(preprocess$network)){
      # browser()
      formula_tmp = update(formula,paste0(formula[[2]],"[[",i,"]]~."))
      terms_to_change = terms[indicator_data]
      if(length(indicator_data)>0){
        for(j in 1:length(terms_to_change)){
          tmp_change = terms_to_change[j]
          name_term = sub(".*=.", "", tmp_change)
          name_term = sub(")", "", name_term)
          name_term_new = paste0(name_term, "[[",i,"]]")
          formula_tmp = update(formula_tmp,new = formula(paste0("~ . - ",
                                                                paste0(c(tmp_change,sub(x = tmp_change,pattern = name_term,replacement = name_term_new)),collapse = " + "))))
        }
      }
      # browser()
      formulae[[i]] = formula_tmp
      preprocess_list[[i]] = formula_preprocess(formula_tmp)
      global_stats[[i]] = count_statistics(formula_tmp)
      
      # loglk = tsergm_path_sampling(preprocess = preprocess, global_stats = global_stats,
      #                              coef = final_coef,coef_indep = final_coef_sub, 
      #                              cluster = cluster,
      #                              briges = control$n_bridges,sampler = control$sampler_path_sampling)
      # # loglk = sergm_path_sampling(formula_tmp, coef = final_coef,coef_indep = final_coef_sub,
      # #                             cluster = control$cluster_bridge,
      # #                             briges = control$n_bridges,sampler = control$sampler_path_sampling)
      # loglk_ls[[i]] = loglk
      # cat(",",i+1)
    }
    if(length(control$cluster) == 0){
      loglk_ls = lapply(X = 1:length(preprocess$network), FUN = function(x){
        tsergm_path_sampling(preprocess = preprocess_list[[x]],global_stats = global_stats[[x]],
                             coef = final_coef,coef_indep = final_coef_sub, briges = control$n_bridges,
                             sampler =control$sampler_path_sampling)
      })
    } else{
      loglk_ls = parLapply(cl = control$cluster, X = 1:length(preprocess$network), fun = function(x){
        tsergm_path_sampling(preprocess = preprocess_list[[x]],global_stats = global_stats[[x]],
                             coef = final_coef,coef_indep = final_coef_sub, briges = control$n_bridges,
                             sampler =control$sampler_path_sampling)
      })
    }
    # 
    # for(i in 1:length(preprocess$network)){
    #   formula_tmp = update(formula,paste0(formula[[2]],"[[",i,"]]~."))
    #   terms_to_change = terms[indicator_data]
    #   if(length(indicator_data)>0){
    #     for(j in 1:length(terms_to_change)){
    #       tmp_change = terms_to_change[j]
    #       name_term = sub(".*=.", "", tmp_change)
    #       name_term = sub(")", "", name_term)
    #       name_term_new = paste0(name_term, "[[",i,"]]")
    #       formula_tmp = update(formula_tmp,new = formula(paste0("~ . - ",
    #                                                             paste0(c(tmp_change,sub(x = tmp_change,pattern = name_term,replacement = name_term_new)),collapse = " + "))))
    #     }
    #   }
    #   # browser()
    #   formulae[[i]] = formula_tmp
    #   preprocess = formula_preprocess(formula_tmp)
    #   global_stats = count_statistics(formula_tmp)
    #   
    #   loglk = tsergm_path_sampling(preprocess = preprocess, global_stats = global_stats,
    #                                coef = final_coef,coef_indep = final_coef_sub, 
    #                                cluster = control$cluster,
    #                               briges = control$n_bridges,sampler = control$sampler_path_sampling)
    #   # loglk = sergm_path_sampling(formula_tmp, coef = final_coef,coef_indep = final_coef_sub,
    #   #                             cluster = control$cluster_bridge,
    #   #                             briges = control$n_bridges,sampler = control$sampler_path_sampling)
    #   loglk_ls[[i]] = loglk
    #   cat(",",i+1)
    # }
    # cat("\n")
    loglk = sum(unlist(loglk_ls))
    # res$loglik = loglk - llk.dind + llk.null
    # browser()
    res$loglik = loglk + sum(log(sub_model$fitted.values^(sub_model$y)))
    
    # res$loglik = loglk + sum(log(sub_model$fitted.values[sub_model$y == 1]))
  }
  res$formula = formula
  res$n_actors = n_actors
  cat("Finished\n")
  class(res) = "tsergm"
  return(res)
}
# export
eval_loglik = function(res, formula, cluster = NULL){
  # browser()
  control = res$control
  # Turn the formula into terms and network data 
  preprocess = t_formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network[[1]])
  final_coef = res$coefficients
  cat("Evaluate the log likelihood\n")
  data_list_dyad = list()
  for(i in 1:length(preprocess$network)){
    data_list_dyad[[i]] = preprocess$data_list[[i]][preprocess$dyad_idep]
  }
  
  sub_model = t_mple_est(networks =preprocess$network,n_actors = n_actors,mod_return = T,
                         terms = preprocess$term_names[preprocess$dyad_idep],
                         data_lists = data_list_dyad,
                         type_list = preprocess$type_list[preprocess$dyad_idep], 
                         coef_names = preprocess$coef_names[preprocess$dyad_idep])
  final_coef_sub = final_coef
  final_coef_sub[] = 0
  # llk.dind = -sub_model$deviance/2 - -sub_model$null.deviance/2
  # llk.null = -log( (n_actors*(n_actors-1)/2)^3)
  final_coef_sub[match( names(sub_model$coefficients),names(final_coef_sub))] = sub_model$coefficients
  loglk_ls = list()
  terms = attr(terms.formula(formula),"term.labels")
  indicator_data = grep(pattern = "data =",x = terms)
  formulae = preprocess_list = global_stats = list()
  for(i in 1:length(preprocess$network)){
    # browser()
    formula_tmp = update(formula,paste0(formula[[2]],"[[",i,"]]~."))
    terms_to_change = terms[indicator_data]
    if(length(indicator_data)>0){
      for(j in 1:length(terms_to_change)){
        tmp_change = terms_to_change[j]
        name_term = sub(".*=.", "", tmp_change)
        name_term = sub(")", "", name_term)
        name_term_new = paste0(name_term, "[[",i,"]]")
        formula_tmp = update(formula_tmp,new = formula(paste0("~ . - ",
                                                              paste0(c(tmp_change,sub(x = tmp_change,pattern = name_term,replacement = name_term_new)),collapse = " + "))))
      }
    }
    # browser()
    formulae[[i]] = formula_tmp
    preprocess_list[[i]] = formula_preprocess(formula_tmp)
    global_stats[[i]] = count_statistics(formula_tmp)
    
    # loglk = tsergm_path_sampling(preprocess = preprocess, global_stats = global_stats,
    #                              coef = final_coef,coef_indep = final_coef_sub, 
    #                              cluster = cluster,
    #                              briges = control$n_bridges,sampler = control$sampler_path_sampling)
    # # loglk = sergm_path_sampling(formula_tmp, coef = final_coef,coef_indep = final_coef_sub,
    # #                             cluster = control$cluster_bridge,
    # #                             briges = control$n_bridges,sampler = control$sampler_path_sampling)
    # loglk_ls[[i]] = loglk
    # cat(",",i+1)
  }
  if(length(cluster) == 0){
    loglk_ls = lapply(X = 1:length(preprocess$network), FUN = function(x){
      tsergm_path_sampling(preprocess = preprocess_list[[x]],global_stats = global_stats[[x]],
                           coef = final_coef,coef_indep = final_coef_sub, briges = control$n_bridges,
                           sampler =control$sampler_path_sampling)
    })
  } else{
    loglk_ls = parLapply(cl = cluster, X = 1:length(preprocess$network), fun = function(x){
      tsergm_path_sampling(preprocess = preprocess_list[[x]],global_stats = global_stats[[x]],
                           coef = final_coef,coef_indep = final_coef_sub, briges = control$n_bridges,
                           sampler =control$sampler_path_sampling)
    })
  }
  
  
  cat("\n")
  loglk = sum(unlist(loglk_ls))
  # res$loglik = loglk - llk.dind + llk.null
  res$loglik = loglk + sum(log(sub_model$fitted.values[sub_model$y == 1]))
  return(res)
}
#' Summary Method for tsergm Objects
#'
#' Prints basic information of the results of an estimated sergm model. 
#'
#' @param object A fitted sergm model (object of class \code{sergm}).
#' @param ... Optional parameter (currently not in use).
#' @return table of the results as a data.frame object.
#' @export
summary.tsergm = function (object, ...) 
{    control = object$control

if(control$method_var == "Fisher"){
  coef = object$coefficients
  ans = list(formula = object$formula, call = object$call, var = object$var, 
             estimate = object$coefficients, 
             control = object$control)
  mod.se = sqrt(diag(object$var))
  est.se = sqrt(diag(object$mcmc_var))
  tot.se = sqrt(diag(object$var+object$mcmc_var))
  est.pct = est.se/tot.se
  zval = coef/tot.se
  pval = 2 * pnorm(q = abs(zval), lower.tail = FALSE)
  count = 1
  coefmat = cbind(Estimate = coef, `Std. Error` = tot.se, `MCMC %` = est.pct,
                  `z value` = zval, `Pr(>|z|)` = pval)
  rownames(coefmat) <- names(object$coefficients)
  if(control$eval_likelihood == T){
    mle.lik = object$loglik
  }
  
  ans$coefficients <- coefmat
} else if (control$method_var == "MPLE"){
  control = object$control
  coef = object$coefficients
  ans = list(formula = object$formula, call = object$call, var = object$var, 
             estimate = object$coefficients, 
             control = object$control)
  mod.se = sqrt(diag(object$var))
  zval = coef/mod.se
  pval = 2 * pnorm(q = abs(zval), lower.tail = FALSE)
  count = 1
  coefmat = cbind(Estimate = coef, `Std. Error` = mod.se,
                  `z value` = zval, `Pr(>|z|)` = pval)
  rownames(coefmat) <- names(object$coefficients)
  if(control$eval_likelihood == T){
    mle.lik = object$loglik
  }
  
  ans$coefficients <- coefmat
}

class(ans) <- "summary.tsergm"
ans
}


#' AIC for tsergm Objects
#'
#' Prints the AIC of an estimated sergm model. 
#'
#' @param object A fitted tsergm model (object of class \code{tsergm}).
#' @param ... optionally more fitted model objects (not used at the moment).
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @return numeric value of AIC.
#' @export
AIC.tsergm = function(object, ..., k = 2) {
  loglik = object$loglik
  df = length(object$coefficients)
  return(k*df-2*loglik)
}

#' Draw from the distribution of an Signed Exponential Family Random Graph Model
#' 
#' \code{sergm_simulate} is used to draw from exponential
#' family random network models for signed networks.  See \code{\link{sergm}} for more
#' information on these models. 
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
#' @param cluster A PSOCK or Fork cluster generated with \code{\link{makeCluster}} from the package \code{\link{parallel}}.
#' 
#' @return Depending on the parameter \code{only_stats}, either a list including the sampled networks 
#' and a matrix of their sufficient statistics or a list only including the matrix of sufficient statistics is returned.  
#' Each network is represented as a list of two lists. One list for the negative edges and one for the positive edges. 
#' This list can be transformed back to a matrix with the function \code{\link{map_to_mat}}.   
#' 
#' @export
tsergm_simulate = function(formula, coef,only_stats =T ,sampler = NULL, cluster = NULL){
  # Set default values 
  if(is.null(sampler)){
    sampler = sampler.sergm()
  }
  # Turn the formula into terms and network data 
  preprocess = t_formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network[[1]])
  if(length(cluster) == 0){
    res = lapply(X = 1:length(preprocess$network), FUN = function(x){
      if(sampler$init_empty == T){
        res = simulate_networks_new(coefs = coef ,n_actors = n_actors,
                                    n_proposals_burn_in = sampler$n_proposals_burn_in,
                                    n_proposals = sampler$n_proposals,seed = sampler$seed,
                                    number_networks = sampler$number_networks,
                                    terms = preprocess$term_names,
                                    data = preprocess$data_list[[x]],
                                    type =preprocess$type_list,only_stats = only_stats,mh = sampler$mh)
        colnames(res$stats) = preprocess$coef_names
      } else {
        res = simulation_mat(network = preprocess$network[[x]], coef = coef ,n_actors = n_actors,
                             n_proposals = sampler$n_proposals,seed = sampler$seed,
                             number_networks = sampler$number_networks,
                             terms = preprocess$term_names,
                             data_list = preprocess$data_list[[x]],
                             type_list =preprocess$type_list,
                             mh = sampler$mh)
        res = rbind(t(res[[1]]),res[[2]])
        colnames(res) = preprocess$coef_names
        res= list(stats = res)
      }
      return(res)
    })
  } else{
    res = parLapply(cl = cluster, X = 1:length(preprocess$network), fun = function(x){
      if(sampler$init_empty == T){
        res = simulate_networks_new(coefs = coef ,n_actors = n_actors,
                                    n_proposals_burn_in = sampler$n_proposals_burn_in,
                                    n_proposals = sampler$n_proposals,seed = sampler$seed,
                                    number_networks = sampler$number_networks,
                                    terms = preprocess$term_names,
                                    data = preprocess$data_list[[x]],
                                    type =preprocess$type_list,only_stats = only_stats,mh = sampler$mh)
        colnames(res$stats) = preprocess$coef_names
      } else {
        res = simulation_mat(network = preprocess$network[[x]], coef = coef ,n_actors = n_actors,
                             n_proposals = sampler$n_proposals,seed = sampler$seed,
                             number_networks = sampler$number_networks,
                             terms = preprocess$term_names,
                             data_list = preprocess$data_list[[x]],
                             type_list =preprocess$type_list,
                             mh = sampler$mh)
        res = rbind(t(res[[1]]),res[[2]])
        colnames(res) = preprocess$coef_names
        res= list(stats = res)
      }
      return(res)
    })
  }
  
  return(res)
}


#' Model assessment for tsergm Objects
#'
#' Samples from the converged object and then calculate the distribution of 
#' edgewise shared friends and enemies as well as positive and negative degree distributions. 
#'
#' @param object A fitted sergm model (object of class \code{\link{sergm}}).
#' @param sampler A fitted sergm model (object of class \code{\link{sampler.sergm}}).
#' @param cluster A PSOCK or Fork cluster generated with \code{\link{makeCluster}} from the package \code{\link{parallel}}.
#' @return numeric value of AIC.
#' @export
t_model_assessment = function(object, sampler = ergm.sign::sampler.sergm(), cluster = NULL){
  if(!"tsergm" %in% class(object)){
    stop("Object has to be of class 'tsergm'!")
  }
  preprocess = t_formula_preprocess(object$formula)
  n_actors = object$n_actors
  # Step 1: Sample networks according to the model 
  sampled_networks = tsergm_simulate(object$formula,coef =object$coefficients,
                                     only_stats = F,sampler = sampler,cluster = cluster)
  
  # Step 2: Transform the sampled networks to matrices
  if(length(cluster) == 0){
    matrix_networks = lapply(sampled_networks,FUN = function(y){
      lapply(y$networks, FUN = function(x){
        map_to_mat(x$edges_pos,x$edges_neg,n_actors)
      })
    })
  } else {
    matrix_networks = parLapply(cl = cluster,sampled_networks,fun = function(y){
      lapply(y$networks, FUN = function(x){
        map_to_mat(x$edges_pos,x$edges_neg,n_actors)
      })
    })
  }
  # Step 2: Get the statistics of the observed network
  count_ese_list = list()
  count_esp_list = list()
  degree_pos_list = list()
  degree_neg_list = list()
  # igraph_data_list = list()
  for(i in 1:length(preprocess$network)){
    number_pos_edges = sum(preprocess$network[[i]] == 1)
    count_ese = count_edgewise_shared_enemies(network = preprocess$network[[i]],n_actors =object$n_actors,n_edges =  number_pos_edges)
    max_ese = max(count_ese)
    count_ese = factor(count_ese, levels = 0:max_ese)
    count_esp = count_edgewise_shared_partner_pos(network = preprocess$network[[i]],n_actors =object$n_actors,n_edges =  number_pos_edges)
    max_esp = max(count_esp)
    count_esp = factor(count_esp, levels = 0:max_esp)
    degree_pos = rowSums(preprocess$network[[i]] == 1)
    max_dp = max(degree_pos)
    degree_pos = factor(degree_pos, levels = 0:max_dp)
    degree_neg = rowSums(preprocess$network[[i]] == -1)
    max_dn = max(degree_neg)
    degree_neg = factor(degree_neg, levels = 0:max_dn)
    # igraph_tmp = mat_to_sign_igraph(preprocess$network[[i]])
    # igraph_data_obs = c("balance_score" = balance_score(igraph_tmp),count_signed_triangles(igraph_tmp))
    
    count_ese_list[[i]] =table(count_ese)
    count_esp_list[[i]] = table(count_esp)
    degree_pos_list[[i]] = table(degree_pos)
    degree_neg_list[[i]] = table(degree_neg)
    # igraph_data_list[[i]] = igraph_data_obs
  }
  
  # Step 3: Get statistics per sample
  if(length(cluster) == 0){
    statistics = lapply(sampled_networks, function(y){
      lapply(y$networks,FUN = function(x){
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
    })
  } else {
    statistics = parLapply(cl = cluster, X = sampled_networks, fun = function(y){
      lapply(y$networks,FUN = function(x){
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
    })
  }
  
  # browser()
  df_ese_list = list()
  df_esp_list = list()
  df_dp_list = list()
  df_dn_list = list()
  # df_igraph_list = list()
  for(i in 1:length(preprocess$network)){
    df_ese = do.call("rbind", lapply(statistics[[i]], FUN = function(x){data.frame(t(as.numeric(x$ese)))}))
    names(df_ese) = 1:ncol(df_ese)-1
    df_esp = do.call("rbind", lapply(statistics[[i]], FUN = function(x){data.frame(t(as.numeric(x$esp)))}))
    names(df_esp) = 1:ncol(df_esp)-1
    df_dp = do.call("rbind", lapply(statistics[[i]], FUN = function(x){data.frame(t(as.numeric(x$dp)))}))
    names(df_dp) = 1:ncol(df_dp)-1
    df_dn = do.call("rbind", lapply(statistics[[i]], FUN = function(x){data.frame(t(as.numeric(x$dn)))}))
    names(df_dn) = 1:ncol(df_dn) -1
    # df_igraph = do.call("rbind", lapply(statistics[[i]], FUN = function(x){data.frame(t(as.numeric(x$igraph_data)))}))
    # names(df_igraph) = names(statistics[[1]]$igraph_data)
    df_ese_list[[i]] = df_ese
    df_esp_list[[i]] = df_esp
    df_dp_list[[i]] = df_dp
    df_dn_list[[i]] = df_dn
    # df_igraph_list[[i]] = df_igraph
  }
  return(list(df_ese_list = df_ese_list,
              df_esp_list = df_esp_list,
              df_dp_list = df_dp_list,
              df_dn_list= df_dn_list,
              # df_igraph_list = df_igraph_list,
              # igraph_data_obs_list = igraph_data_list,
              count_ese = count_ese_list,
              count_esp = count_esp_list,
              count_dp = degree_pos_list,
              count_dn = degree_neg_list,
              sampler = sampler))
}

