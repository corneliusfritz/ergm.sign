#' @importFrom parallel parLapply
#' @importFrom stats binomial glm qnorm quantile rnorm terms terms.formula
#' @importFrom lpSolveAPI make.lp set.column set.objfn set.constr.type set.rhs set.bounds get.objective lp.control
#' @importFrom stringr str_split
#' @importFrom coda mcmc thin as.mcmc.list mcmc.list
#' @importFrom stats acf ar cov na.pass pnorm update
#' @importFrom Rcpp compileAttributes 
#' @importFrom RcppArmadillo fastLm
# Taken from ergm see https://github.com/statnet/ergm/blob/fd83fc94bf077a962d753966fe53d51627cb3136/R/is.inCH.R
is.inCH <- function(p, M, verbose=FALSE, ...) { # Pass extra arguments directly to LP solver
  p = t(p)
  verbose <- max(0, min(verbose, 4))
  
  if(is.null(dim(p))) p <- rbind(p)
  
  if (!is.matrix(M)) 
    stop("Second argument must be a matrix.")
  if ((d <- ncol(p)) != ncol(M))
    stop("Number of columns in matrix (2nd argument) is not equal to dimension ",
         "of first argument.")
  
  if((n <- nrow(M)) == 1L){
    for(i in seq_len(nrow(p))){
      if(!isTRUE(all.equal(p[i,], M, check.attributes = FALSE))) return(FALSE)
    }
    return(TRUE)
  }
  
  # Set up the optimisation problem: the following are common for all rows of p.
  
  timeout <- 1
  
  setup.lp <- function(){
    L <- cbind(1, M)
    lprec <- make.lp(n, d+1)
    for(k in seq_len(d+1)) set.column(lprec, k, L[,k])
    set.constr.type(lprec, rep.int(2L, n)) # 2 = ">="
    set.rhs(lprec,  numeric(n))
    set.bounds(lprec, lower = rep.int(-1, d+1L), upper = rep.int(1, d+1L))
    lp.control(lprec, break.at.value = -.Machine$double.eps, verbose=c("important","important","important","normal","detailed")[min(max(verbose+1,0),5)], timeout=timeout)
    lprec
  }
  lprec <- setup.lp()
  
  for(i in seq_len(nrow(p))){ # Iterate over test points.
    
    # Keep trying until results are satisfactory.
    #
    # flag meanings:
    # -1      : dummy value, just starting out
    #  0 or 11: Good (either 0 or some negative value)
    #  1 or  7: Timeout
    #   others: probably nothing good, but don't know how to handle
    flag <- -1
    while(flag%in%c(-1,1,7)){
      # Set the objective function in terms of p and solve the problem.
      set.objfn(lprec, c(1, p[i,]))
      flag <- solve(lprec)
      if(flag%in%c(1,7)){ # Timeout
        timeout <- timeout * 2 # Increase timeout, in case it's just a big problem.
        z <- rnorm(1) # Shift target and test set by the same constant.
        p <- p + z
        M <- M + z
        lprec <- setup.lp() # Reinitialize
      }
    }
    
    # If the objective function (min) is not zero, the point p[i,] is not in the CH of M.
    if(get.objective(lprec) < 0){
      if(verbose > 1) message(sprintf("is.inCH: iter = %d, outside hull.",i))
      return(FALSE)
    }else if(verbose > 2 && nrow(p) > 1) message(sprintf("is.inCH: iter = %d, inside hull.", i))
    
  }
  
  if(verbose > 1) message("is.inCH: all test points inside hull.")
  return(TRUE) # If all points passed the test, return TRUE.
}
#' Function to find the gamma in the stepping algorithm as high as possible while still being in the convex hull of the sampled statistics 
#' @param gammas Vector of gamma values to be tried out 
#' @param obs_statistics Matrix of simulated statistics (columns the separate sufficient statistics and rows the separate simulations)
#' @param M_bar Average simulated statistics over all simulations 
#' @param M Observed simulated statistics

#' @export
find_min_in_ch = function(gammas, obs_statistics, M, M_bar){
  tmp = FALSE
  i = length(gammas)+1
  obs_statistics = t(obs_statistics)
  # M_bar = t(M_bar)
  while(tmp == FALSE){
    i = i -1
    tmp = is.inCH(p = as.vector(1.05*gammas[i]*obs_statistics + (1-1.05*gammas[i])*M_bar),M = M)
  }
  return(i-1)
  
}
# export
formula_preprocess = function(formula = network ~ edges_neg + edges_pos){
  network = eval(formula[[2]])
  # network = as.matrix(as_adjacency_matrix(network,sparse = F,attr = "sign"))
  
  tmp_type_vec = c()
  relating_vec = c()
  tmp_type_names = c()
  tmp_type_vec_only_data = c()
  tmp_type_names_only_data = c()
  
  terms = attr(terms.formula(formula),"term.labels")
  # if(length(terms) == 0){
  #   terms = as.character(formula[[3]])
  # }
  # Where special terms including data or type arguments used? 
  indicator_type = grep(pattern = "type =",x = terms)
  indicator_data = grep(pattern = "data =",x = terms)
  indicator_only_data = indicator_data[which(!indicator_data %in% indicator_type)]
  indicator_normal = which(!(1:length(terms)) %in% c(indicator_type, indicator_data))
  data_list_of_typed_terms = list()
  data_list_of_only_data_terms = list()
  n = 1
  if(length(indicator_type)>0){
    terms_with_type = terms[indicator_type]
    for(i in 1:length(terms_with_type)) {
      name_term = sub("type.*", "", terms_with_type[i]) 
      name_term = gsub(pattern = "\\(",replacement = "",x = name_term)
      tmp_term = str_split(terms_with_type[i]," ")[[1]]
      # Where is type in cut_term?
      cut_term = tmp_term[grep(pattern = "type",tmp_term) + 2]
      # Does the term also include data? 
      data_term = tmp_term[grep(pattern = "data",tmp_term) + 2]
      data_term = gsub(pattern = ",",replacement = " ",x = data_term)
      data_term = trimws(gsub(pattern = ")",replacement = " ",x = data_term))
      # cut_term = sub(".*type =", "", terms_with_type[i]) 
      cut_term = gsub(pattern = "))",replacement = ")",x = cut_term)
      cut_term = gsub(pattern = ",",replacement = " ",x = cut_term)
      tmp_type = eval(parse(text = cut_term))
      names(tmp_type) = paste(name_term,data_term,tmp_type,sep = "_")
      tmp_type_vec = c(tmp_type_vec,tmp_type)
      tmp_type_names = c(tmp_type_names,rep(name_term, length(tmp_type)))
      
      relating_vec = c(relating_vec,   
                       match(terms_with_type[i],terms) + seq(from = 0, to = 0.5,length.out = length(tmp_type)))
      
      if(length(data_term) >0){
        for(j in 1:length(tmp_type)){
          data_list_of_typed_terms[[n]] =  eval(parse(text = data_term))
          n = n +1
        }
      } else {
        data_list_of_typed_terms = c(data_list_of_typed_terms, replicate(length(tmp_type_names), matrix(data = 0), simplify=FALSE))
        n = n + length(tmp_type_names)
      }
    }
  }
  
  n = 1
  tmp = 0
  if(length(indicator_only_data)>0){
    terms_only_data = terms[indicator_only_data]
    for(i in 1:length(terms_only_data)) {
      name_term = sub("data.*", "", terms_only_data[i]) 
      name_term = gsub(pattern = "\\(",replacement = "",x = name_term)
      tmp_term = str_split(terms_only_data[i]," ")[[1]]
      data_term = tmp_term[grep(pattern = "data",tmp_term) + 2]
      data_term = trimws(gsub(pattern = ")",replacement = " ",x = data_term))
      # Add data
      data_list_of_only_data_terms[[n]] =  eval(parse(text = data_term))
      names(tmp) = paste(name_term, data_term, sep = "_")
      tmp_type_vec_only_data = c(tmp_type_vec_only_data,tmp)
      tmp_type_names_only_data = c(tmp_type_names_only_data,name_term)
      relating_vec = c(relating_vec,   
                       match(terms_only_data[i],terms) )
      
      n = n + 1
    }
  }
  
  data_per_term = c(data_list_of_typed_terms,data_list_of_only_data_terms,replicate(length(indicator_normal), matrix(data = 0), simplify=FALSE))
  type_per_term = c(tmp_type_vec, tmp_type_vec_only_data, vector("numeric", length =  length(indicator_normal)))
  coef_name_per_term = c(names(tmp_type_vec), names(tmp_type_vec_only_data),terms[indicator_normal])
  coef_name_per_term = sub(pattern = "\\[\\[",replacement = "_",x = coef_name_per_term)
  coef_name_per_term = sub(pattern = "\\]\\]",replacement = "",x = coef_name_per_term)
  name_per_term = c(tmp_type_names,tmp_type_names_only_data,terms[indicator_normal])
  # relating_vec = c(relating_vec, terms[indicator_normal])
  relating_vec = c(relating_vec,   
                   match(terms[indicator_normal],terms) )
  coef_name_per_term = gsub(pattern = "__",replacement = "_",x = coef_name_per_term)
  dyad_idep = name_per_term %in% c("edges_neg", "edges_pos",
                                   "edges", "cov_dyad_pos", 
                                   "cov_dyad_neg", "cov_dyad")
  return(list(network = network, terms = terms,
              data_list = data_per_term[order(relating_vec)], 
              type_list = type_per_term[order(relating_vec)], 
              coef_names = coef_name_per_term[order(relating_vec)], 
              term_names = name_per_term[order(relating_vec)], 
              dyad_idep = dyad_idep[order(relating_vec)]))
}
# export
t_formula_preprocess = function(formula = networks ~ edges_neg + edges_pos + cov_dyad_pos(data = cov_pos)+ cov_dyad_neg(data = cov_pos)){
  networks = eval(formula[[2]])
  terms = attr(terms.formula(formula),"term.labels")
  indicator_data = grep(pattern = "data =",x = terms)
  indicator_no_data = (1:length(terms))[! (1:length(terms)) %in% indicator_data]
  formulae = list()
  
  for(i in 1:length(networks)){
    tmp_formula = formula(formula)
    
    terms_to_change = terms[indicator_data]
    # For each entry with data we change the data argument to the j-th entry of the data argument 
    if(length(indicator_data)>0){
      for(j in 1:length(terms_to_change)){
        tmp_change = terms_to_change[j]
        name_term = sub(".*=.", "", tmp_change)
        name_term = sub(")", "", name_term)
        name_term_new = paste0(name_term, "[[",i,"]]")
        tmp_formula = update(tmp_formula,new = formula(paste0("~ . - ",
                                                              paste0(c(tmp_change,sub(x = tmp_change,pattern = name_term,replacement = name_term_new)),collapse = " + "))))
      }
    }
    formulae[[i]] =  formula_preprocess(tmp_formula)
  }
  data_lists = lapply(formulae, FUN = function(x)x$data_list)
  tmp_formula = formula(formula)
  
  if(length(terms_to_change)>0){
    for(j in 1:length(terms_to_change)){
      tmp_change = terms_to_change[j]
      name_term = sub(".*=.", "", tmp_change)
      name_term = sub(")", "", name_term)
      name_term_new = name_term
      tmp_formula = update(tmp_formula,new = formula(paste0("~ . - ",
                                                            paste0(c(tmp_change,sub(x = tmp_change,pattern = name_term,replacement = name_term_new)),collapse = " + "))))
    }    
  } 
  base_preprocess = formula_preprocess(tmp_formula)
  return(list(network = networks, 
              terms = base_preprocess$terms,
              data_list = data_lists, 
              type_list = base_preprocess$type_list, 
              coef_names = base_preprocess$coef_names, 
              term_names = base_preprocess$term_names, 
              dyad_idep = base_preprocess$dyad_idep))
  
}

#' Sampler for Signed Exponential Random Graph Models 
#'
#' \code{sampler.sergm}  is used to specify MCMC samplers to sample signed networks 
#' from the distribution of Signed Exponential Random Graph Models. These samplers 
#' are needed to carry out MCMC estimation routines, quantify the uncertainty of the  
#' estimates and assess the goodness-of-fit of the estimated models. 
#' 
#' @param n_proposals_burn_in Number of iterations of the sampler to be carried out 
#' as the burn-in period. This number should be picked to be high enough to reach the 
#' stationary distribution of the SERGM from which one wants to sample from.    
#' @param n_proposals Number of iterations of the sampler to be carried out 
#' after the burn-in period for each network (often called MCMC interval).
#' @param seed Seed for the samples to make the code reproducible. 
#' @param number_networks Number of networks to sample. 
#' @param init_empty Boolean value indicating whether the initial network should be the observed network or an empty network.
#' @param mh Boolean value indicating whether an Metropolis Hastings or Gibbs Sampler should be used. 
#' If the value is TRUE the Metropolis Hastings algorithm is used, otherwise the Gibbs Sampler. 
#' @return Object of class \code{sampler.sergm}.
#' @export
sampler.sergm = function(n_proposals_burn_in = 10000*20, n_proposals = 10000,
                         seed = NA,number_networks = 100,init_empty = T, mh = T) {
  if(is.na(seed)){
    seed = sample(1:1000,size = 1)
  }
  res = list(n_proposals_burn_in = n_proposals_burn_in, n_proposals = n_proposals,
             seed = seed,number_networks = number_networks,init_empty = init_empty, mh = mh) 
  class(res) = "sampler.sergm"
  return(res)
}

#' Control the specifications of the estimation algorithm of a Signed Exponential Random Graph Models 
#'
#' \code{control.sergm}  is used to specify how the estimation of a Signed 
#' Exponential Random Graph Models should be carried out.  
#' 
#' @param method_est There are currently four estimation algorithms implemented. \code{MLE} is the standard MCMCLE estimation routine 
#' with a log-normal approximation as proposed by Hummel et al. for binary networks. 
#' \code{Stepping} implements the stepping algorithm also introduced by by Hummel et al., 
#' while \code{RM} evokes the Robins Monroe algorithm for signed networks. 
#' Finally, the pseudo likelihood estimates based on a multinomial approximation 
#' of the intractable likelihood can be obtain through the option \code{MPLE}. 
#' For the pseudo likelihood the standard errors are approximated by a parametric 
#' bootstrap procedure introduced by Cranmer and Desmarais for static binary ERGMs. 
#' @param method_var Method to be used to assess the uncertainty of the estimates. There are three possible options 
#' to be chosen from. The first option is "Fisher" and the standard choice, being based on approximating 
#' the Fisher information via the sampler given in the argument \code{sampler_var}. Second, 
#' "BS" translates to the bootstrap approximation detailed in Desmarais and Cranmer (2012). 
#' Finally, the variance of the MPLE can be used when the parameter equals "MPLE", although here some caution is in place as they are 
#' often too slim and their behavior is not studied in theoretical contexts. If no parameter is specified the default 
#' paramter is the "BS" option for estimations carried out via "MPLE" and "Fisher" for all MCMC-based estimation routines.
#' @param eval_likelihood Boolean indicator whether the logarithmic likelihood should be evaluated or not.
#' @param mple_init Boolean indicator whether the MPLE estimator or a vector of 
#' zeroes should be used as the initial estimator. 
#' @param max_it Maximal number of iterations of the iterative procedure. 
#' @param tol Define the tolerance of relative change in parameter vector between consecutive 
#' iterations to indicate convergence of the estimates. 
#' @param RM_c If method = "RM" the gain sequence a_i in the i-th iteration is given by c^(-i). 
#' @param Stepping_number_grids If method = "Stepping" the number of grids between 0 and 1 for which it will be tested 
#' whether the convex combination between the observed sufficient statistics and the mean simulated statistics lies 
#' in the convex hull of the samples of sufficient statistics.   
#' @param n_bridges Number of bridges used for evaluating the log likelihood via path sampling. 
#' @param sampler_est An \code{\link{sampler.sergm}} to specify the sampler used for the estimation of the parameters 
#' This parameter is not of importance if the chosen method is the MPLE. 
#' @param sampler_var An \code{\link{sampler.sergm}} to specify the sampler used for assessing the uncertainty
#' of the parameters of the estimated \code{sergm}.
#' @param sampler_path_sampling An \code{\link{sampler.sergm}} to specify the sampler used for evaluating the log likelihood at the estimated parameter value. 
#' If the sampler is set to NULL the configurations of \code{sampler_var} are copied. 
#' @param cluster A PSOCK or Fork cluster generated with \code{\link{makeCluster}} from the package \code{\link{parallel}}.
#' @return Object of class \code{sampler.sergm}.
#' @export
control.sergm = function(method_est = "Stepping", method_var = NULL, eval_likelihood = TRUE, 
                         mple_init = TRUE, max_it = 10, tol = 0.003,
                         RM_c = 0.5,Stepping_number_grids = 20, n_bridges = 16, 
                         sampler_est = NULL, 
                         sampler_var = NULL, 
                         sampler_path_sampling = NULL,
                         cluster = NULL) { 
  
  # Set default values 
  if(is.null(sampler_est)){
    sampler_est = ergm.sign::sampler.sergm(init_empty = TRUE, n_proposals_burn_in = 10000, 
                                           n_proposals = 10000,mh = FALSE)
  }
  if(is.null(sampler_var)){
    sampler_var = ergm.sign::sampler.sergm(init_empty = FALSE,
                                           n_proposals_burn_in = 0,
                                           n_proposals = 2000,
                                           number_networks = 1000, 
                                           mh = TRUE)
  }
  if(is.null(sampler_path_sampling)){
    sampler_path_sampling = sampler_var
  }
  # If no seed was given pick a random integer between 1 and 1000 (do it for the var and est sampler)
  if(is.null(sampler_est$seed)){
    sampler_est$seed = sample(1:1000,size = 1)
  }
  if(is.null(sampler_var$seed)){
    sampler_var$seed = sample(1:1000,size = 1)
  }
  if(!method_est %in% c("Stepping", "MLE", "RM", "MPLE")){
    stop("The estimator needs to be one of the following: Stepping, MLE, RM, MPLE")
  }
  # Set the default values according to method_est as detailed in the documentation
  if(is.null(method_var)){
    if(method_est == "MPLE"){
      method_var = "BS"
    } else {
      method_var = "Fisher"
    }
  }
  if((tol<= 0) | (tol>1)){
    stop("The parameter tol has to be a numeric value bigger than 0 and smaller or equal than 1.")
  }
  if((RM_c>1) |(RM_c<0.5)){
    stop("The parameter RM_c has to be a numeric value between 0.5 and 1.")
  }
  if((mple_init != TRUE) & (mple_init != FALSE)){
    stop("The parameter mple_init has to be a boolean value.")
  }
  if((eval_likelihood != TRUE) & (eval_likelihood != FALSE)){
    stop("The parameter mple_init has to be a boolean value.")
  }
  if(Stepping_number_grids%%1!=0){
    stop("The parameter Stepping_number_grids has to be a natural number.")
  }
  if(n_bridges%%1!=0){
    stop("The parameter n_bridges has to be a natural number.")
  }
  if(class(sampler_est) != "sampler.sergm"){
    stop("sampler_est needs to be a sampler.sergm object.")
  }
  if(class(sampler_var) != "sampler.sergm"){
    stop("sampler_var needs to be a sampler.sergm object.")
  }
  
  if(class(sampler_path_sampling) != "sampler.sergm"){
    stop("sampler_path_sampling needs to be a sampler.sergm object.")
  }
  if(!method_var %in% c("Fisher", "BS", "MPLE")){
    stop("The variance estimator needs to be one of the following: Fisher, BS, MPLE")
  }
  res = list(method_est = method_est,method_var = method_var, 
             eval_likelihood = eval_likelihood,
             mple_init = mple_init, max_it = max_it, tol = tol, 
             RM_c = RM_c, Stepping_number_grids = Stepping_number_grids, 
             n_bridges = n_bridges, 
             sampler_est = sampler_est, 
             sampler_var = sampler_var, 
             sampler_path_sampling = sampler_path_sampling, 
             cluster = cluster)
  class(res) = "control.sergm"
  return(res)
}
# This function simulates but only gives out all proposals and the changes they evoke 
sergm_simulate_trying = function(formula, coef,only_stats =T ,control = sampler.sergm()){
  # Turn the formula into terms and network data 
  preprocess = formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network)
  # if(length(coef) != preprocess$term_names){
  #   stop("Number of coefficients is not equal the statistics!")
  # }
  res = simulate_trying(coefs = coef ,n_actors = n_actors,
                        n_proposals_burn_in = control$n_proposals_burn_in,
                        n_proposals = control$n_proposals,seed = control$seed,
                        number_networks = control$number_networks,
                        terms = preprocess$term_names,
                        data = preprocess$data_list,
                        type =preprocess$type_list,only_stats = only_stats)
  return(res)
}


#' Convert a matrix to an igraph object 
#'
#' \code{mat_to_sign_igraph}  is used to covert an adjacency matrix representing a signed network to an igraph object.   
#' 
#' @param mat numeric matrix 
#' @return igraph network.
mat_to_sign_igraph = function(mat){
  edges_pos = data.frame(which(mat== 1, arr.ind = T))
  edges_pos$sign = 1
  edges_neg = data.frame(which(mat== -1, arr.ind = T))
  edges_neg$sign = -1
  # edges_neg$sign = factor(edges_neg$sign)
  n_actors = nrow(mat)
  return(tbl_graph(edges = rbind(edges_pos, edges_neg),nodes = data.frame(name = 1:n_actors),directed = F))
}

#' Convert a map to a matrix
#'
#' \code{map_to_mat}  is used to covert a map (which is a c++ data structure) equal to a list to an adjacency matrix representing a signed network.
#' 
#' @param edges_pos List of length n_actors that gives the positive edges for each actor
#' @param edges_neg List of length n_actors that gives the negative edges for each actor
#' @param n_actors Numeric value indicating the number of actors in the network
#' @return igraph network.
map_to_mat = function(edges_pos, edges_neg,n_actors){
  # Generate empty network
  network_tmp = matrix(0,nrow = n_actors,ncol = n_actors)
  for(i in 1:n_actors){
    network_tmp[as.numeric(names(edges_pos)[i]),edges_pos[[i]]] = 1
    network_tmp[as.numeric(names(edges_neg)[i]),edges_neg[[i]]] = -1
  }
  # class(network_tmp) = "snetwork"
  return(network_tmp)
}

#' Calculate the global statistics given in the formula in the signed network on the lhs
#'
#' 
#' @param formula formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a numeric matrix including 
#'   the signed entries 0, 1, and -1. Which terms are implemented is detailed in \code{\link{sergm.terms}}
#' @export
count_statistics = function(formula){
  preprocess = formula_preprocess(formula)
  n_actors = nrow(preprocess$network)
  res = as.vector(count_global_terms(network = preprocess$network,terms = preprocess$term_names,
                                     n_actors = n_actors,data_list = preprocess$data_list, type_list = preprocess$type_list))
  names(res) = preprocess$coef_names
  return(res)
}

#' Calculate the global statistics given in the formula in the signed network on the lhs
#'
#' 
#' @param formula formula An \R \code{\link{formula}} object, of the form
#'   \code{y ~ <model terms>}, where \code{y} is a numeric matrix including 
#'   the signed entries 0, 1, and -1. Which terms are implemented is detailed in \code{\link{sergm.terms}}
#' @export
temporal_count_statistics = function(formula){
  # Turn the formula into terms and network data 
  preprocess = t_formula_preprocess(formula)
  # Get the number of actors 
  n_actors = nrow(preprocess$network[[1]])
  
  res = matrix(nrow = length(preprocess$network),ncol = length(preprocess$terms))
  for(i in 1:length(preprocess$network)){
    res[i,] = as.vector(count_global_terms(network = preprocess$network[[i]],terms = preprocess$term_names,
                                 n_actors = n_actors,data_list = preprocess$data_list[[i]], type_list = preprocess$type_list))
  }
  colnames(res) = preprocess$coef_names
  return(res)
}

get_change_stats = function(formula){
  preprocess = formula_preprocess(formula)
  n_actors = nrow(preprocess$network)
  data_tmp = preprocess_pseudo_lh_new(network = preprocess$network,
                                      n_actors = n_actors,
                                      terms = preprocess$term_names,
                                      data_list = preprocess$data_list,type_list = preprocess$type_list)
  number_terms = length(terms)
  y = c(as.numeric(data_tmp$pos_response), 
        as.numeric(data_tmp$neg_response), 
        as.numeric(data_tmp$zero_response))
  data_processed = data.frame(y)
  for(i in 1:number_terms){
    eval(parse(text = paste0("data_processed$",preprocess$coef_names[i],
                             "= c(data_tmp$pos_changes[[",i,"]], data_tmp$neg_changes[[",i,"]], data_tmp$zero_changes[[",i,"]])")))
  }
  
  # tmp_data = data.table(y,x_1 = x_1-global_stats[1],x_2 = x_2-global_stats[2])
  data_processed = data_processed[!is.nan(data_processed[[2]]),]
  # data_processed = data_processed[,.(weight = .N),by = c(names(data_processed))]
  return(data_processed)
}


NVL = function (...) {
  for (e in eval(substitute(alist(...)))) {
    x <- eval(e, parent.frame())
    if (!is.null(x)) 
      break
  }
  x
}

.catchToList = function (expr) 
{
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), 
                  error = eHandler)
  list(value = val, warnings = myWarnings, error = myError)
}

ERRVL = function (...) 
{
  x <- NULL
  for (e in eval(substitute(alist(...)))) {
    x <- eval(if (inherits(x, "try-error")) 
      do.call(substitute, list(e, list(. = x)))
      else e, parent.frame())
    if (!inherits(x, "try-error")) 
      return(x)
  }
  stop("No non-error expressions passed.")
}

NVL2 = function (test, notnull, null = NULL) 
{
  if (is.null(test)) 
    null
  else notnull
}
# export
sergm.MCMCse <- function (theta0, theta,statsmatrices,H) { 
  av <- colMeans(statsmatrices)
  xsims <- sweep(statsmatrices, 2, av)
  # to change
  gsims <--xsims
  xobs <- -av
  xsim <- as.matrix(xsims)
  gsim <- as.matrix(gsims)
  prob <- rep.int(1/nrow(xsim), nrow(xsim))
  cov.zbar <- spectrum0.mvar(gsims) * sum(prob^2)
  imp.factor <- sum(prob^2) * length(prob)
  novar <- rep(TRUE, nrow(H))
  novar <- diag(H) < sqrt(.Machine$double.eps)
  cov.zbar.obs <- cov.zbar
  cov.zbar.obs[, ] <- 0
  H.obs <- H
  H.obs[, ] <- 0
  imp.factor.obs <- NULL
  H <- H[!novar, !novar, drop = FALSE]
  H.obs <- H.obs[!novar, !novar, drop = FALSE]
  cov.zbar <- cov.zbar[!novar, !novar, drop = FALSE]
  cov.zbar.obs <- cov.zbar.obs[!novar, !novar, drop = FALSE]
  mc.cov.offset <- matrix(0, ncol = length(theta), nrow = length(theta))
  H <- H.obs - H
  cov.zbar <- cov.zbar + cov.zbar.obs
  mc.cov <- matrix(NA, ncol = length(novar), nrow = length(novar))
  mc.cov0 <- solve(H, cov.zbar, tol = 1e-20)
  mc.cov0 <- solve(H, t(mc.cov0), tol = 1e-20)
  mc.cov[!novar, !novar] <- mc.cov0
  attr(mc.cov, "imp.factor") <- imp.factor
  return(mc.cov)
}


spectrum0.mvar = function (x, order.max = NULL, aic = is.null(order.max), tol = .Machine$double.eps^0.5, ...) { 
  breaks <- if (is.mcmc.list(x)) 
    c(0, cumsum(sapply(x, niter)))
  else NULL
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  v <- matrix(0, p, p)
  novar <- abs(apply(x, 2, stats::sd)) < tol
  x <- x[, !novar, drop = FALSE]
  if (ncol(x) == 0) 
    stop("All variables are constant.")
  first_local_min <- function(x) {
    d <- diff(c(Inf, x, Inf))
    min(which(d >= 0)) - 1
  }
  e <- eigen(cov(x), symmetric = TRUE)
  Q <- e$vectors[, sqrt(pmax(e$values, 0)/max(e$values)) > 
                   tol * 2, drop = FALSE]
  xr <- x %*% Q
  ind.var <- cov(xr)
  xr <- if (!is.null(breaks)) 
    do.call(mcmc.list, lapply(lapply(seq_along(breaks[-1]), 
                                     function(i) xr[(breaks[i] + 1):(breaks[i + 1]), , 
                                                    drop = FALSE]), mcmc))
  else as.mcmc.list(mcmc(xr))
  ord <- NVL(order.max, ceiling(10 * log10(niter(xr))))
  xr <- do.call(rbind, c(lapply(unclass(xr)[-nchain(xr)], function(z) rbind(cbind(z), 
                                                                            matrix(NA, ord, nvar(z)))), unclass(xr)[nchain(xr)]))
  arfit <- .catchToList(ar(xr, aic = is.null(order.max), order.max = ord, 
                           na.action = na.pass, ...))
  while ((!is.null(arfit$error) || ERRVL(try(any(eigen(arfit$value$var.pred, 
                                                       only.values = TRUE)$values < 0), silent = TRUE), TRUE)) && 
         ord > 0) {
    ord <- ord - 1
    if (ord <= 0) 
      stop("Unable to fit ar() even with order 1; this is likely to be due to insufficient sample size or a trend in the data.")
    arfit <- .catchToList(ar(xr, aic = is.null(order.max), 
                             order.max = ord, na.action = na.pass, ...))
  }
  arfit <- arfit$value
  if (aic && arfit$order > (ord <- first_local_min(arfit$aic) - 
                            1)) {
    arfit <- ar(xr, aic = ord == 0, order.max = max(ord, 
                                                    1), na.action = na.pass)
  }
  arvar <- arfit$var.pred
  arcoefs <- arfit$ar
  arcoefs <- NVL2(dim(arcoefs), apply(arcoefs, 2:3, base::sum), 
                  sum(arcoefs))
  adj <- diag(1, nrow = ncol(xr)) - arcoefs
  iadj <- solve(adj)
  v.var <- iadj %*% arvar %*% t(iadj)
  infl <- exp((determinant(v.var)$modulus - determinant(ind.var)$modulus)/ncol(ind.var))
  v.var <- Q %*% v.var %*% t(Q)
  v[!novar, !novar] <- v.var
  attr(v, "infl") <- infl
  v
}
