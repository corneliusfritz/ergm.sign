#' Parallel Sampling 
#'
#' \code{\link{parallel_sample}}  is used to sample temporal networks in parallel. This is an internal function and not to be used directly.  
#'
#' @param tmp_coef numeric vector of coefficients to be sampled from 
#' @param terms set of terms 
#' @param n_actors Numeric value of the number of actors
#' @param networks List of signed networks (each represeting a snaptshot in time) 
#' @param data_lists List of data matrices for each term
#' @param type_list List of type arguments for each term
#' @param mh boolean value indicating whether Metropolis-Hastings or Gibbs sampling should be used
#' @param n_proposals Number of iterations between sampled networks (number of tries)
#' @param n_proposals_burn_in Number of burn in iterations
#' @param seed Numeric seed to make the stochastic sampling reproducible 
#' @param number_networks Numeric value of the number of networks to be sampled. 
#' @param cluster A PSOCK or Fork cluster generated with \code{\link{makeCluster}} from the package \code{\link{parallel}}.
#' @param start_with_empty_net boolean value indicating whether we should start with the empty or observed network
#' @return Object of class \code{sergm}.
#' @export
parallel_sample= function(tmp_coef,terms, n_actors,
                          networks, data_lists, type_list,mh,
                          n_proposals,n_proposals_burn_in,
                          seed,number_networks, cluster,start_with_empty_net) {
  if(start_with_empty_net){
    sampled_stats = parLapply(cl = cluster, X = 1:length(data_lists),fun = function(x){
      simulate_networks_stats(tmp_coef,terms, n_actors,
                              data_lists[[x]], type_list,mh,
                              n_proposals,n_proposals_burn_in,
                              seed + x,number_networks)
    })
  } else {
    sampled_stats = parLapply(cl = cluster, X = 1:length(data_lists),fun = function(x){
      simulation_mat(networks[[x]],terms, n_actors,n_proposals,seed + x,number_networks,
                     data_lists[[x]], type_list,tmp_coef,mh)[[2]]
    })
  }
  
  sampled_means = lapply(sampled_stats,FUN = function(x){colMeans(x)})
  
  sampled_vars = lapply(sampled_stats,FUN = function(x){cov(x)})
  sum_means = Reduce("+",sampled_means)
  sum_vars = Reduce("+",sampled_vars)
  sum_stats = Reduce("+",sampled_stats)
  
  return(list(sum_means,sum_vars,sum_stats))
}

mple_est = function(network,n_actors = n_actors, mod_return = F,terms,data_list,type_list, coef_names) {
  data_tmp = preprocess_pseudo_lh_new(network,n_actors = n_actors,terms = terms,data_list = data_list,type_list = type_list)
  number_terms = length(terms)
  y = c(as.numeric(data_tmp$pos_response), 
        as.numeric(data_tmp$neg_response), 
        as.numeric(data_tmp$zero_response))
  data_processed = data.frame(y)
  for(i in 1:number_terms){
    eval(parse(text = paste0("data_processed$",coef_names[i],
                             "= c(data_tmp$pos_changes[[",i,"]], data_tmp$neg_changes[[",i,"]], data_tmp$zero_changes[[",i,"]])")))
  }
  
  # tmp_data = data.table(y,x_1 = x_1-global_stats[1],x_2 = x_2-global_stats[2])
  data_processed = data_processed[!is.nan(data_processed[[2]]),]
  # data_processed = data_processed[,.(weight = .N),by = c(names(data_processed))]
  tmp_model = glm(y~ .-1, family = binomial(), data = data_processed)
  if(mod_return){
    return(tmp_model)
  } else {
    return(as.numeric(tmp_model$coefficients))
  }
}

t_mple_est = function(networks,n_actors = n_actors, mod_return = F,terms,data_lists,type_list, coef_names) {
  data_processed_list = list()
  for(j in 1:length(networks)){
    data_tmp = preprocess_pseudo_lh_new(networks[[j]],n_actors = n_actors,terms = terms,data_list = data_lists[[j]],
                                        type_list = type_list)
    number_terms = length(terms)
    y = c(as.numeric(data_tmp$pos_response), 
          as.numeric(data_tmp$neg_response), 
          as.numeric(data_tmp$zero_response))
    data_processed = data.frame(y)
    for(i in 1:number_terms){
      eval(parse(text = paste0("data_processed$",coef_names[i],
                               "= c(data_tmp$pos_changes[[",i,"]], data_tmp$neg_changes[[",i,"]], data_tmp$zero_changes[[",i,"]])")))
    }
    # tmp_data = data.table(y,x_1 = x_1-global_stats[1],x_2 = x_2-global_stats[2])
    data_processed_list[[j]] = data_processed[!is.nan(data_processed[[2]]),]
  }
  data_processed = do.call("rbind", data_processed_list)
  tmp_model = glm(y~ .-1, family = binomial(), data = data_processed)
  if(mod_return){
    return(tmp_model)
  } else {
    return(as.numeric(tmp_model$coefficients))
  }
}
# Function that carries out the path sampling
# The function bears some similarities to the ERGM implementation 
# export
sergm_path_sampling = function(formula, coef,coef_indep, briges = 20,sampler = sampler.sergm(),cluster = NA){
  # Taken from the ergm package 
  from = coef_indep
  to = coef
  # Taken from ergm see https://github.com/statnet/ergm/blob/master/R/ergm.bridge.R
  mkpath <- function(n, shift = 0, reverse = FALSE) {
    stopifnot(shift >= -1/2, shift <= 1/2)
    u0 <- seq(from = 0 + 1 / 2 / n, to = 1 - 1 / 2 / n, length.out = n)
    u <- u0 + shift / n
    if (reverse) u <- rev(u)
    list(
      theta = t(rbind(sapply(u, function(u) cbind(to * u + from * (1 - u))))),
      u = u
    )
  }
  
  path <- mkpath(briges, shift = 0, reverse = F)
  Dtheta.Du <- (to-from)/briges
  # Turn the formula into terms and network data 
  preprocess = formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network)
  global_stats = count_statistics(formula)
  llrs = numeric(length = briges)
  # cat("Using ", briges, " bridges: ")
  if(length(cluster) == 0){
    llrs = lapply(X = seq_len(briges), FUN = function(x){
      theta <- path$theta[x, ]
      if(sampler$init_empty == T){
        res = simulate_networks_stats(coefs = theta ,n_actors = n_actors,
                                      n_proposals_burn_in = sampler$n_proposals_burn_in,
                                      n_proposals = sampler$n_proposals,seed = sampler$seed+x,
                                      number_networks = sampler$number_networks,
                                      terms = preprocess$term_names,
                                      data = preprocess$data_list,
                                      type =preprocess$type_list,mh = sampler$mh)
      } else {
        res = simulation_mat(network = preprocess$network, coef = coef ,n_actors = n_actors,
                             n_proposals = sampler$n_proposals,seed = sampler$seed,
                             number_networks = sampler$number_networks,
                             terms = preprocess$term_names,
                             data_list = preprocess$data_list,
                             type_list =preprocess$type_list,
                             mh = sampler$mh)
        res = rbind(t(res[[1]]),res[[2]])
      }
      
      return(mean(res %*%Dtheta.Du))
      # return(mean(-res %*%Dtheta.Du) - (-global_stats%*% Dtheta.Du))
    })
  } else{
    llrs = parLapply(cl = cluster, X = seq_len(briges), fun = function(x){
      theta <- path$theta[x, ]
      if(sampler$init_empty == T){
        res = simulate_networks_stats(coefs = theta ,n_actors = n_actors,
                                      n_proposals_burn_in = sampler$n_proposals_burn_in,
                                      n_proposals = sampler$n_proposals,seed = sampler$seed+x,
                                      number_networks = sampler$number_networks,
                                      terms = preprocess$term_names,
                                      data = preprocess$data_list,
                                      type =preprocess$type_list,mh = sampler$mh)
      } else {
        res = simulation_mat(network = preprocess$network, coef = coef ,n_actors = n_actors,
                             n_proposals = sampler$n_proposals,seed = sampler$seed,
                             number_networks = sampler$number_networks,
                             terms = preprocess$term_names,
                             data_list = preprocess$data_list,
                             type_list =preprocess$type_list,
                             mh = sampler$mh)
        res = rbind(t(res[[1]]),res[[2]])
      }
      return(mean(res %*%Dtheta.Du))
      # return(mean(-res %*%Dtheta.Du) - (-global_stats%*% Dtheta.Du))
    })
  }
  return((to-from)%*%global_stats - sum(unlist(llrs)))
  
  # return(sum(unlist(llrs)))
}

# export


tsergm_path_sampling = function(preprocess, global_stats,coef,coef_indep, briges = 20,sampler = sampler.sergm(),cluster = NULL){
  # Taken from the ergm package 
  from = coef_indep
  to = coef
  # Taken from ergm see https://github.com/statnet/ergm/blob/master/R/ergm.bridge.R
  mkpath <- function(n, shift = 0, reverse = FALSE) {
    stopifnot(shift >= -1/2, shift <= 1/2)
    u0 <- seq(from = 0 + 1 / 2 / n, to = 1 - 1 / 2 / n, length.out = n)
    u <- u0 + shift / n
    if (reverse) u <- rev(u)
    list(
      theta = t(rbind(sapply(u, function(u) cbind(to * u + from * (1 - u))))),
      u = u
    )
  }
  
  path <- mkpath(briges, shift = 0, reverse = F)
  Dtheta.Du <- (to-from)/briges
  # Turn the formula into terms and network data 
  # Get the number of actors 
  n_actors = nrow(preprocess$network)
  llrs = numeric(length = briges)
  # cat("Using ", briges, " bridges: ")
  
  tmp_fun = function(x,path,n_actors,sampler,preprocess) {
    theta <- path$theta[x, ]
    if(sampler$init_empty == T){
      res = simulate_networks_stats(coefs = theta ,n_actors = n_actors,
                                    n_proposals_burn_in = sampler$n_proposals_burn_in,
                                    n_proposals = sampler$n_proposals,seed = sampler$seed+x,
                                    number_networks = sampler$number_networks,
                                    terms = preprocess$term_names,
                                    data = preprocess$data_list,
                                    type =preprocess$type_list,mh = sampler$mh)
    } else {
      res = simulation_mat(network = preprocess$network, coef = coef ,n_actors = n_actors,
                           n_proposals = sampler$n_proposals,seed = sampler$seed,
                           number_networks = sampler$number_networks,
                           terms = preprocess$term_names,
                           data_list = preprocess$data_list,
                           type_list =preprocess$type_list,
                           mh = sampler$mh)
      res = rbind(t(res[[1]]),res[[2]])
    }
    
    return(mean(res %*%Dtheta.Du))
  }
  if(length(cluster) == 0){
    llrs = lapply(X = seq_len(briges), FUN = tmp_fun,path = path,n_actors = n_actors,
                  sampler = sampler,preprocess = preprocess)
  } else{
    llrs = parLapply(cl = cluster, X = seq_len(briges), fun = tmp_fun,path = path,n_actors = n_actors,
                     sampler = sampler,preprocess = preprocess)
  }
  return((to-from)%*%global_stats - sum(unlist(llrs)))
  
  # return(sum(unlist(llrs)))
}