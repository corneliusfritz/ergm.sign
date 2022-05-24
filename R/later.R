
# Obtain tie probability for non-edge-exchangeable models
#
# Samples from the converged object and then calculates the relative frequency of each dyadwise outcome.
#
# @param object A fitted sergm model (object of class \code{\link{sergm}}).
# @param sampler A fitted sergm model (object of class \code{\link{sampler.sergm}}).
# @return numeric value of AIC.
# export
get_tie_probability = function(object, sampler = ergm.sign::sampler.sergm()){
  if(!"sergm" %in% class(object)){
    stop("Object has to be of class 'sergm'!")
  }
  preprocess = formula_preprocess(object$formula)
  n_actors = nrow(preprocess$network)
  # Step 1: Sample networks according to the model 
  sampled_networks = sergm_simulate(object$formula,coef =object$coefficients,only_stats = F,sampler = sampler)
  # Step 2: Transform the sampled networks to matrices
  matrix_networks = lapply(sampled_networks$networks, FUN = function(x){
    map_to_mat(x$edges_pos,x$edges_neg,n_actors)
  })
  matrix_networks_pos = lapply(sampled_networks$networks, FUN = function(x){
    map_to_mat(x$edges_pos,x$edges_neg,n_actors) == 1
  })
  matrix_networks_neg = lapply(sampled_networks$networks, FUN = function(x){
    map_to_mat(x$edges_pos,x$edges_neg,n_actors) == -1
  })
  perc_pos = 1/sampler$number_networks*Reduce("+",matrix_networks_pos)
  perc_neg = 1/sampler$number_networks*Reduce("+",matrix_networks_neg)
  perc_zero = 1 - perc_pos - perc_neg
  
  pred = matrix(data = NA,nrow = n_actors,ncol = n_actors)
  pred[(perc_pos>perc_neg) &(perc_pos>perc_zero)] = 1
  pred[(perc_neg>perc_pos) &(perc_neg>perc_zero)] = -1
  pred[(perc_zero>perc_neg) &(perc_zero>perc_pos)] = 0
  diag(preprocess$network) = NA
  
  acc = sum(pred == preprocess$network, na.rm = T)/(n_actors*(n_actors-1))
  
  confusion_matrix_sample = matrix(data = c(  sum((preprocess$network == 1)*perc_pos, na.rm = T),
                                              sum((preprocess$network == 1)*perc_neg, na.rm = T),
                                              sum((preprocess$network == 1)*perc_zero, na.rm = T),
                                              sum((preprocess$network == -1)*perc_pos, na.rm = T),
                                              sum((preprocess$network == -1)*perc_neg, na.rm = T),
                                              sum((preprocess$network == -1)*perc_zero, na.rm = T),
                                              sum((preprocess$network == 0)*perc_pos, na.rm = T),
                                              sum((preprocess$network == 0)*perc_neg, na.rm = T),
                                              sum((preprocess$network == 0)*perc_zero, na.rm = T)),nrow = 3,ncol = 3)
  
  
  corr_per_truth = confusion_matrix_sample/colSums(confusion_matrix_sample)
  corr_per_pred = confusion_matrix_sample/rowSums(confusion_matrix_sample)
  return(list(acc = acc, confusion_matrix = confusion_matrix_sample, 
              corr_per_truth = corr_per_truth, 
              corr_per_pred = corr_per_pred))
  
}


# Obtain tie probability for non-edge-exchangeable models
#
# Samples from the converged object and then calculates the relative frequency of each dyadwise outcome.
#
# @param object A fitted sergm model (object of class \code{\link{sergm}}).
# @param cluster A PSOCK or Fork cluster generated with \code{\link{makeCluster}} from the package \code{\link{parallel}}.
# @param sampler A fitted sergm model (object of class \code{\link{sampler.sergm}}).
# @return numeric value of AIC.
# export
t_get_tie_probability = function(object, cluster = NULL, 
                                 sampler = ergm.sign::sampler.sergm()){
  if(!"tsergm" %in% class(object)){
    stop("Object has to be of class 'sergm'!")
  }
  preprocess = t_formula_preprocess(object$formula)
  n_actors = nrow(preprocess$network[[1]])
  # Step 1: Sample networks according to the model 
  sampled_networks = tsergm_simulate(object$formula,coef =object$coefficients,
                                     only_stats = F,sampler = sampler)
  
  # Step 2: Transform the sampled networks to matrices
  matrix_networks = lapply(sampled_networks,FUN = function(y){
    lapply(y$networks, FUN = function(x){
      map_to_mat(x$edges_pos,x$edges_neg,n_actors)
    })
  })
  
  matrix_networks_pos = lapply(sampled_networks,FUN = function(y){
    lapply(y$networks, FUN = function(x){
      map_to_mat(x$edges_pos,x$edges_neg,n_actors) == 1
    })
  })
  matrix_networks_neg =lapply(sampled_networks,FUN = function(y){
    lapply(y$networks, FUN = function(x){
      map_to_mat(x$edges_pos,x$edges_neg,n_actors) == -1
    })
  })
  
  perc_pos =lapply(matrix_networks_pos, FUN = function(x){
    1/sampler$number_networks*Reduce("+",x)
  })
  perc_neg =lapply(matrix_networks_neg, FUN = function(x){
    1/sampler$number_networks*Reduce("+",x)
  })
  perc_zero =lapply(1:length(matrix_networks_neg), FUN = function(x){
    1 - perc_pos[[x]] - perc_neg[[x]]
  })
  pred_list = list()
  confusion_matrix_sample = list()
  corr_per_truth = list()
  corr_per_pred = list()
  acc = c()
  for(i in 1:length(matrix_networks_neg)){
    pred = matrix(data = NA,nrow = n_actors,ncol = n_actors)
    pred[(perc_pos[[i]]>perc_neg[[i]]) &(perc_pos[[i]]>perc_zero[[i]])] = 1
    pred[(perc_neg[[i]]>perc_pos[[i]]) &(perc_neg[[i]]>perc_zero[[i]])] = -1
    pred[(perc_zero[[i]]>perc_neg[[i]]) &(perc_zero[[i]]>perc_pos[[i]])] = 0
    diag(pred) = NA
    acc[[i]] = sum(pred == preprocess$network[[i]], na.rm = T)/(n_actors*(n_actors-1))
    
    confusion_matrix_sample[[i]] = matrix(data = c(  sum((preprocess$network[[i]] == 1)*perc_pos[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == 1)*perc_neg[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == 1)*perc_zero[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == -1)*perc_pos[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == -1)*perc_neg[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == -1)*perc_zero[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == 0)*perc_pos[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == 0)*perc_neg[[i]], na.rm = T),
                                                     sum((preprocess$network[[i]] == 0)*perc_zero[[i]], na.rm = T)),nrow = 3,ncol = 3)
    corr_per_truth[[i]] = confusion_matrix_sample[[i]]/colSums(confusion_matrix_sample[[i]])
    corr_per_pred[[i]] = confusion_matrix_sample[[i]]/rowSums(confusion_matrix_sample[[i]])
    
  }
  
  
  
  
  
  
  return(list(acc = acc, confusion_matrix = confusion_matrix_sample, 
              corr_per_truth = corr_per_truth, 
              corr_per_pred = corr_per_pred, 
              perc_pos = perc_pos, 
              perc_neg = perc_neg, 
              perc_zero = perc_zero))
}
