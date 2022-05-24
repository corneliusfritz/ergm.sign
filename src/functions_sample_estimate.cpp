// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <random>
#include <unordered_set>
#include <unordered_map>    
#include <set>
#include <map> 
#include <list> 
#include "change_statistics.h"
#include "function_helper.h"
#include <RcppArmadilloExtensions/sample.h>
// # include <omp.h> 
#include <unistd.h>
//// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]] 
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec arma_sample(int n, int k, arma::vec prob){
  arma::uvec out(n, arma::fill::zeros);
  arma::vec comp = arma::regspace(1,  k);
  // arma::vec comp = {1,-1,0}; 
  prob.print();
  for(int i = 0; i < n; ++i){
    out(i) = Rcpp::RcppArmadillo::sample(comp, 1, false, prob)[0];
  }
  prob.print();
  return out;
}
// [[Rcpp::export]]
int sample_signed(int n, arma::vec prob){
  arma::vec comp = {1,-1,0};
  return Rcpp::RcppArmadillo::sample(comp, 1, false, prob)[0];
}

arma::vec eval_at_empty_network(std::vector<std::string> terms, int n_actors){
  arma::vec res(terms.size());
  for(unsigned int i = 0; i < terms.size(); i += 1){
    // Add here any additional change statistics 
    if(terms.at(i) == "edges_pos") {
      res.at(i) = 0;
    } else if(terms.at(i) == "edges_neg") {
      res.at(i) = 0;
    } else if(terms.at(i) == "cov_dyad_pos") {
      res.at(i) = 0;
    } else if(terms.at(i) == "cov_dyad_neg") {  
      res.at(i) = 0;
    } else if(terms.at(i) == "two_star_pos") {
      res.at(i) = 0;
    } else if(terms.at(i) == "two_star_neg") {
      res.at(i) = 0;
    } else if(terms.at(i) == "degree_pos") { 
      res.at(i) = 0;
    } else if(terms.at(i) == "degree_neg") {
      res.at(i) = 0;
    } else if(terms.at(i) == "gwdegree_pos") { 
      res.at(i) = 0;  
    } else if(terms.at(i) == "gwdegree_neg") { 
      res.at(i) = 0;
    } else if(terms.at(i) == "gwese_pos") { 
      res.at(i)= 0;    
    } else if(terms.at(i) == "gwese_neg") { 
      res.at(i)= 0;    
    } else if(terms.at(i) == "gwdsp_pos") { 
      res.at(i)= 0;    
    } else if(terms.at(i) == "gwesf_pos") { 
      res.at(i)= 0;     
    } else if(terms.at(i) == "gwesf_neg") { 
      res.at(i)= 0;     
    }else if(terms.at(i) == "gwdasp") { 
      res.at(i)= 0;     
    } else if(terms.at(i) == "gwdegree") { 
      res.at(i)= 0;     
    } else if(terms.at(i) == "degree") { 
      res.at(i)= 0;      
    } else if(terms.at(i) == "edges") { 
      res.at(i)= 0;     
    } else if(terms.at(i) == "cov_dyad") { 
      res.at(i)= 0;     
    } else if(terms.at(i) == "gwdsp_neg") { 
      res.at(i)= 0;     
    } else if(terms.at(i) == "isolates"){
      res.at(i)= n_actors;     
    }else if(terms.at(i) == "isolates_pos"){
      res.at(i)= n_actors;     
    }else if(terms.at(i) == "isolates_neg"){
      res.at(i)= n_actors;     
    }else if(terms.at(i) == "cov_dyad_neut"){
      res.at(i)= 0;     
    }
  }
  return(res);
}

// Function to generate vector of functions according to the vector of 
std::vector<arma::vec(*)(const  std::unordered_map< int, std::set<int>> &edges_pos,
                      const std::unordered_map< int, std::set<int>> &edges_neg,
                      int &n_actors,
                      int &tmp_random_i,
                      int &tmp_random_j, 
                      arma::mat &data, 
                      int &type)> change_statistics_generate(std::vector<std::string> terms){
                        std::vector<arma::vec(*)( const  std::unordered_map< int, std::set<int>> &edges_pos,
                                              const  std::unordered_map< int, std::set<int>> &edges_neg,
                                              int &n_actors,
                                              int &tmp_random_i,
                                              int &tmp_random_j, 
                                              arma::mat &data, 
                                              int &type)> functions;
                        // Now we need to find for each term the right formula 
                        for(unsigned int i = 0; i < terms.size(); i += 1){
                          // Add here any additional change statistics 
                          if(terms.at(i) == "edges_pos") {
                            functions.push_back(stat_edges_pos);
                          } else if(terms.at(i) == "edges_neg") {
                            functions.push_back(stat_edges_neg);
                          } else if(terms.at(i) == "cov_dyad_pos") {
                            functions.push_back(stat_cov_dyad_pos);
                          } else if(terms.at(i) == "cov_dyad_neg") {  
                            functions.push_back(stat_cov_dyad_neg);
                          } else if(terms.at(i) == "two_star_pos") {
                            functions.push_back(stat_two_star_pos);  
                          } else if(terms.at(i) == "two_star_neg") {
                            functions.push_back(stat_two_star_neg); 
                          } else if(terms.at(i) == "degree_pos") {
                            functions.push_back(stat_degree_pos);
                          } else if(terms.at(i) == "degree_neg") {
                            functions.push_back(stat_degree_neg);
                          } else if(terms.at(i) == "gwdegree_pos") {   
                            functions.push_back(stat_gwdegree_pos);    
                          } else if(terms.at(i) == "gwdegree_neg") {  
                            functions.push_back(stat_gwdegree_neg);
                          } else if(terms.at(i) == "gwesf_pos") { 
                            functions.push_back(stat_gwesf_pos);
                          } else if(terms.at(i) == "gwesf_neg") { 
                            functions.push_back(stat_gwesf_neg);
                          } else if(terms.at(i) == "gwdsp_pos") { 
                            functions.push_back(stat_gwdsp_pos);
                          } else if(terms.at(i) == "gwese_pos") {        
                            functions.push_back(stat_gwese_pos);
                          } else if(terms.at(i) == "gwese_neg") {        
                            functions.push_back(stat_gwese_neg);
                          } else if(terms.at(i) == "gwdasp") {        
                            functions.push_back(stat_gwdasp);
                          } else if(terms.at(i) == "gwdegree") {        
                            functions.push_back(stat_gwdegree);
                          } else if(terms.at(i) == "degree") {        
                            functions.push_back(stat_degree);
                          } else if(terms.at(i) == "edges") {        
                            functions.push_back(stat_edges);
                          } else if(terms.at(i) == "cov_dyad") {        
                            functions.push_back(stat_cov_dyad);
                          } else if(terms.at(i) == "gwdsp_neg") { 
                            functions.push_back(stat_gwdsp_neg);
                          } else if(terms.at(i) == "isolates_pos") { 
                            functions.push_back(stat_isolates_pos);
                          }else if(terms.at(i) == "isolates_neg") { 
                            functions.push_back(stat_isolates_neg);
                          }else if(terms.at(i) == "isolates") { 
                            functions.push_back(stat_isolates);
                          }else if(terms.at(i) == "cov_dyad_neut"){
                            functions.push_back(stat_cov_dyad_neut);
                          }
                          
                        }
                        return(functions);
                      }
// Function to call all functions in the vector functions 
arma::mat calculate_change_stats( int &tmp_random_i,
                                  int tmp_random_j, 
                                  int n_actors,
                                  std::unordered_map< int, std::set<int>> &edges_pos,
                                  std::unordered_map< int, std::set<int>> &edges_neg,
                                  std::vector<arma::mat> &data_list, 
                                  std::vector<int> &type_list,
                                  std::vector<arma::vec(*)( const  std::unordered_map< int, std::set<int>> &edges_pos,
                                                        const  std::unordered_map< int, std::set<int>> &edges_neg,
                                                        int &n_actors,
                                                        int &tmp_random_i,
                                                        int &tmp_random_j, 
                                                        arma::mat &data, 
                                                        int &type)> functions){
  arma::mat change_stat(functions.size(),3);
  arma::vec res_tmp;
  for(unsigned int a = 0; a <(functions.size()); a = a + 1 ) {
    // res_tmp = functions.at(a)(n_actors,present_val,proposed_change);
    // Rcout << functions.size() << std::endl;
    // Rcout << data_list.size() << std::endl;
    // Rcout << type_list.size() << std::endl;
    
    change_stat.row(a) = functions.at(a)(edges_pos, edges_neg, n_actors,
                    tmp_random_i,tmp_random_j, data_list.at(a), type_list.at(a)).t();
    // Rcout << "Changed from 1 to 0" << std::endl;
  }
  return(change_stat);
}


// // [[Rcpp::export]]
// arma::mat calculate( int &proposed_change,
//                                     int &present_val, int n_actors,
//                                    std::vector<std::string> terms){
//   arma::mat change_stat(terms.size(),3);
//   // std::list<arma::vec(*)(int,int,int)> functions;
//   std::vector<arma::vec(*)( std::unordered_map< int, std::unordered_set<int>> &edges_pos,
//                          std::unordered_map< int, std::unordered_set<int>> &edges_neg,
//                          int n_actors,
//                          int present_val,
//                          int proposed_change)> functions;
//   functions = change_statistics_generate(terms);
//   arma::vec res_tmp;
//   for(int a = 0; a <(functions.size()); a = a + 1 ) {
//     // res_tmp = functions.at(a)(n_actors,present_val,proposed_change);
//     change_stat.row(a) = functions.at(a)(edges_pos, edges_neg, n_actors,present_val,proposed_change).t();
//   }
//   return(change_stat);
// }

// This is a function that picks the right row of the change stats for the proposed change 
// (this is needed since the change statistics are calculating for each term three change stats 0 to 1, 0 to -1, and 1 to -1)
arma::vec pick_row( int present_val,
                    int proposed_change) {
  arma::vec res(2); // {sign,row}
  if(present_val == 1){
    if(proposed_change == 0){
      res = {-1,1};
    } else {
      res = {1,3};
    }
  }
  if(present_val == -1){
    if(proposed_change == 0){
      res = {-1,2};
    } else {
      res = {-1,3};
    }
  }
  if(present_val == 0){
    if(proposed_change == 1){
      res = {1,1};
    } else {
      res = {1,2};
    }
  }
  return(res);
}

int helper_fun(int proposed_change) {
  if(proposed_change == 1){
    return(0);
  } else if(proposed_change == -1){
    return(1);
  } else {
    return(2);
  } 
}

// Write a function that takes the global_stats, present_val (of edge),
// matrix of change stats to derive (delta(y_ij^+), delta(y_ij^-), delta(y_ij^0))
arma::mat get_global_stats(int& present_val, arma::mat& change_stat) {
  arma::mat res(change_stat.n_rows,3); 
  arma::vec zeros(change_stat.n_rows); 
  zeros.fill(0);
  if(present_val == 1){
    res.col(0) = zeros;
    res.col(1) =  change_stat.col(2);
    res.col(2) = - change_stat.col(0);
  }
  if(present_val == -1){     
    res.col(0) = - change_stat.col(2);
    res.col(1) =  zeros;
    res.col(2) = - change_stat.col(1);
  }
  if(present_val == 0){
    res.col(0) = change_stat.col(0);
    res.col(1) = change_stat.col(1);
    res.col(2) =zeros;
  }
  return(res);
}

// One sweep for simulating a network through n_proposals is carried out via Gibbs Sampling
void simulate_network_gs_tnt( arma::vec coefs,
                              std::unordered_map< int, std::set<int>> &edges_pos,
                              std::unordered_map< int, std::set<int>> &edges_neg,
                              int &n_actors, 
                              int &n_proposals, 
                              int seed, 
                              std::vector<arma::mat> &data_list, 
                              std::vector<int> &type_list,
                              std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                                                    const std::unordered_map< int, std::set<int>> &edges_neg,
                                                    int &n_actors,
                                                    int &tmp_random_i,
                                                    int &tmp_random_j, 
                                                    arma::mat &data, 
                                                    int &type)> &functions, 
                                                    arma::vec &global_stats) {
  // Set up objects 
  // Plan is to implement TNT sampling in the future ... 
  // where edges are selected with prob. 0.5
  // std::set<char> nonzero_edges;
  int present_val, proposed_change;
  // double MR;
  arma::mat HR;
  // change_staCt = arma::vec(n_elem);
  arma::mat change_stat;
  set_seed(seed);
  arma::vec tmp_stat;
  arma::umat comp_vec;
  // Rcout << n_proposals_effective << std::endl;
  
  // Rcpp::Clock clock;
  int tmp_random_i,tmp_random_j; 
  arma::vec tmp_row, prob;
  arma::vec comp = {1,-1,0};  
  arma::mat tmp;
  arma::vec zeros(functions.size()); 
  zeros.fill(0);
  NumericVector random_accept= runif(n_proposals,0,1);
  std::set< std::pair<int, int>> realized_edges;
  // Go through a loop for each proposed change
  for(int a = 0; a <=(n_proposals-1); a = a + 1 ) {
    // Get random proposal
    if(random_accept.at(a)>0.5){
      
    }
    // tmp_random_i = random_i.at(a);
    // tmp_random_j = random_j.at(a);
    // Whats the present val?
    if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
      present_val = 1;
    } else if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
      present_val = -1;
    } else {
      present_val = 0;
    }
    // clock.tock("Get proposal");
    // 2. Step: Calculate the change statistics 
    // Rcout << "Starting change stats" << std::endl;
    // clock.tick("Calculate change stats");
    change_stat = calculate_change_stats(tmp_random_i,
                                         tmp_random_j,
                                         n_actors,
                                         edges_pos,
                                         edges_neg,  
                                         data_list, 
                                         type_list,
                                         functions);
    // clock.tock("Calculate change stats");
    // clock.tick("Calculate probabilities");
    tmp = get_global_stats(present_val,change_stat);
    prob = {exp(coefs.t()*tmp.col(0)).at(0,0), exp(coefs.t()*tmp.col(1)).at(0,0), exp(coefs.t()*tmp.col(2)).at(0,0)};
    prob = prob/(sum(prob));
    // Rcout << prob << std::endl;
    
    // clock.tock("Calculate probabilities");
    
    // Rcout << tmp << std::endl;
    // clock.tick("Sample");
    proposed_change = Rcpp::RcppArmadillo::sample(comp, 1, false, prob)[0];
    // Rcout << proposed_change << std::endl;
    
    // clock.tock("Sample");
    // clock.tick("Update change");
    // if(present_val != proposed_change){
    //   tmp_row = pick_row(present_val,proposed_change);
    //   tmp_stat=change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0);
    // } else {
    //   tmp_stat = zeros;
    // }
    // Rcout << "New iteration" << std::endl;
    
    // Rcout << present_val << std::endl;
    // Rcout << proposed_change << std::endl;
    // Rcout << tmp.col(helper_fun(proposed_change)) << std::endl;
    global_stats += tmp.col(helper_fun(proposed_change));
    // Rcout << "Start" << std::endl;
    // Rcout << change_stat << std::endl;
    // Rcout << "From " + std::to_string(present_val) + " to " + std::to_string(proposed_change) << std::endl;
    // Rcout << "Between " + std::to_string(tmp_random_i) + " and " + std::to_string(tmp_random_j) << std::endl;     
    // Rcout << tmp.col(helper_fun(proposed_change)) << std::endl;
    // Rcout << global_stats << std::endl;
    // arma::uvec ids = find(sign(global_stats) == -1);
    // if(ids.size()>0){
    //   Rcout << "Start" << std::endl;
    //   Rcout << change_stat << std::endl;
    //   Rcout << "From " + std::to_string(present_val) + " to " + std::to_string(proposed_change) << std::endl;
    //   Rcout << "Between " + std::to_string(tmp_random_i) + " and " + std::to_string(tmp_random_j) << std::endl;
    //   Rcout << tmp.col(helper_fun(proposed_change)) << std::endl;
    //   Rcout << global_stats << std::endl;
    //   break;
    // }
    // global_stats += tmp_stat;
    // clock.tock("Update change");
    // clock.tick("Update");
    
    // Here we modify the network
    if(proposed_change == 0){
      if(present_val == 1){
        // Rcout << "Changed from 1 to 0" << std::endl;
        edges_pos.at(tmp_random_i).erase(tmp_random_j);
        edges_pos.at(tmp_random_j).erase(tmp_random_i);
      } else {
        // Rcout << "Changed from -1 to 0" << std::endl;
        edges_neg.at(tmp_random_i).erase(tmp_random_j);
        edges_neg.at(tmp_random_j).erase(tmp_random_i);
      }
    }
    if(proposed_change == 1){
      if(present_val == -1){
        // Rcout << "Changed from -1 to 1" << std::endl;
        edges_neg.at(tmp_random_i).erase(tmp_random_j);
        edges_neg.at(tmp_random_j).erase(tmp_random_i);
        edges_pos.at(tmp_random_i).insert(tmp_random_j);
        edges_pos.at(tmp_random_j).insert(tmp_random_i);
      } else {
        // Rcout << "Changed from -1 to 0" << std::endl;
        edges_pos.at(tmp_random_i).insert(tmp_random_j);
        edges_pos.at(tmp_random_j).insert(tmp_random_i);
      }
    }
    if(proposed_change == -1){
      if(present_val == 1){
        // Rcout << "Changed from 1 to -1" << std::endl;
        edges_pos.at(tmp_random_i).erase(tmp_random_j);
        edges_pos.at(tmp_random_j).erase(tmp_random_i);
        edges_neg.at(tmp_random_i).insert(tmp_random_j);
        edges_neg.at(tmp_random_j).insert(tmp_random_i);
      } else {
        // Rcout << "Changed from 0 to -1" << std::endl;
        edges_neg.at(tmp_random_i).insert(tmp_random_j);
        edges_neg.at(tmp_random_j).insert(tmp_random_i);
      }
    }
    // clock.tock("Update");
    
    // clock.tock("Change network");
    
  }
  // clock.stop("clock");
}

// One sweep for simulating a network through n_proposals is carried out via Gibbs Sampling
void simulate_network_gs( arma::vec coefs,
                          std::unordered_map< int, std::set<int>> &edges_pos,
                          std::unordered_map< int, std::set<int>> &edges_neg,
                          int &n_actors, 
                          int &n_proposals, 
                          int seed, 
                          std::vector<arma::mat> &data_list, 
                          std::vector<int> &type_list,
                          std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                                                const std::unordered_map< int, std::set<int>> &edges_neg,
                                                int &n_actors,
                                                int &tmp_random_i,
                                                int &tmp_random_j, 
                                                arma::mat &data, 
                                                int &type)> &functions, 
                                                arma::vec &global_stats) {
  // Set up objects 
  // Plan is to implement TNT sampling in the future ... 
  // where edges are selected with prob. 0.5
  // std::set<char> nonzero_edges;
  int present_val, proposed_change;
  // double MR;
  arma::mat HR;
  // change_staCt = arma::vec(n_elem);
  arma::mat change_stat;
  set_seed(seed);
  arma::vec tmp_vec, random_i(n_proposals), random_j(n_proposals), tmp_stat;
  arma::umat comp_vec;
  random_i = simulate_numbers(0,n_actors-1, n_proposals)+1;
  random_j = simulate_numbers(0,n_actors-1, n_proposals)+1;
  // Where is random_i bigger than random_j? 
  comp_vec = (random_i>random_j);
  arma::uvec ids = find(comp_vec == 1); // Find indices
  tmp_vec = random_i(ids);
  random_i(ids) = random_j(ids);
  random_j(ids) = tmp_vec;
  comp_vec = (random_i==random_j);
  ids = find(comp_vec == 0); // Find indices
  random_i = random_i(ids);
  random_j = random_j(ids);
  int n_proposals_effective = random_j.size();
  // Rcout << n_proposals_effective << std::endl;
  
  // Rcpp::Clock clock;
  int tmp_random_i,tmp_random_j; 
  arma::vec tmp_row, prob;
  arma::vec comp = {1,-1,0};  
  arma::mat tmp;
  arma::vec zeros(functions.size()); 
  zeros.fill(0);
  
  // Go through a loop for each proposed change
  for(int a = 0; a <=(n_proposals_effective-1); a = a + 1 ) {
    // Get random proposal
    // clock.tick("Get proposal");
    tmp_random_i = random_i.at(a);
    tmp_random_j = random_j.at(a);
    // Whats the present val?
    if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
      present_val = 1;
    } else if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
      present_val = -1;
    } else {
      present_val = 0;
    }
    // clock.tock("Get proposal");
    // 2. Step: Calculate the change statistics 
    // Rcout << "Starting change stats" << std::endl;
    // clock.tick("Calculate change stats");
    change_stat = calculate_change_stats(tmp_random_i,
                                         tmp_random_j,
                                         n_actors,
                                         edges_pos,
                                         edges_neg,  
                                         data_list, 
                                         type_list,
                                         functions);
    // clock.tock("Calculate change stats");
    // clock.tick("Calculate probabilities");
    tmp = get_global_stats(present_val,change_stat);
    prob = {exp(coefs.t()*tmp.col(0)).at(0,0), exp(coefs.t()*tmp.col(1)).at(0,0), exp(coefs.t()*tmp.col(2)).at(0,0)};
    prob = prob/(sum(prob));
    // Rcout << prob << std::endl;
    
    // clock.tock("Calculate probabilities");
    
    // Rcout << tmp << std::endl;
    // clock.tick("Sample");
    proposed_change = Rcpp::RcppArmadillo::sample(comp, 1, false, prob)[0];
    // Rcout << proposed_change << std::endl;
    
    // clock.tock("Sample");
    // clock.tick("Update change");
    // if(present_val != proposed_change){
    //   tmp_row = pick_row(present_val,proposed_change);
    //   tmp_stat=change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0);
    // } else {
    //   tmp_stat = zeros;
    // }
    // Rcout << "New iteration" << std::endl;
    
    // Rcout << present_val << std::endl;
    // Rcout << proposed_change << std::endl;
    // Rcout << tmp.col(helper_fun(proposed_change)) << std::endl;
    global_stats += tmp.col(helper_fun(proposed_change));
    // Rcout << "Start" << std::endl;
    // Rcout << change_stat << std::endl;
    // Rcout << "From " + std::to_string(present_val) + " to " + std::to_string(proposed_change) << std::endl;
    // Rcout << "Between " + std::to_string(tmp_random_i) + " and " + std::to_string(tmp_random_j) << std::endl;     
    // Rcout << tmp.col(helper_fun(proposed_change)) << std::endl;
    // Rcout << global_stats << std::endl;
    // arma::uvec ids = find(sign(global_stats) == -1);
    // if(ids.size()>0){
    //   Rcout << "Start" << std::endl;
    //   Rcout << change_stat << std::endl;
    //   Rcout << "From " + std::to_string(present_val) + " to " + std::to_string(proposed_change) << std::endl;
    //   Rcout << "Between " + std::to_string(tmp_random_i) + " and " + std::to_string(tmp_random_j) << std::endl;
    //   Rcout << tmp.col(helper_fun(proposed_change)) << std::endl;
    //   Rcout << global_stats << std::endl;
    //   break;
    // }
    // global_stats += tmp_stat;
    // clock.tock("Update change");
    // clock.tick("Update");
    
    // Here we modify the network
    if(proposed_change == 0){
      if(present_val == 1){
        // Rcout << "Changed from 1 to 0" << std::endl;
        edges_pos.at(tmp_random_i).erase(tmp_random_j);
        edges_pos.at(tmp_random_j).erase(tmp_random_i);
      } else {
        // Rcout << "Changed from -1 to 0" << std::endl;
        edges_neg.at(tmp_random_i).erase(tmp_random_j);
        edges_neg.at(tmp_random_j).erase(tmp_random_i);
      }
    }
    if(proposed_change == 1){
      if(present_val == -1){
        // Rcout << "Changed from -1 to 1" << std::endl;
        edges_neg.at(tmp_random_i).erase(tmp_random_j);
        edges_neg.at(tmp_random_j).erase(tmp_random_i);
        edges_pos.at(tmp_random_i).insert(tmp_random_j);
        edges_pos.at(tmp_random_j).insert(tmp_random_i);
      } else {
        // Rcout << "Changed from -1 to 0" << std::endl;
        edges_pos.at(tmp_random_i).insert(tmp_random_j);
        edges_pos.at(tmp_random_j).insert(tmp_random_i);
      }
    }
    if(proposed_change == -1){
      if(present_val == 1){
        // Rcout << "Changed from 1 to -1" << std::endl;
        edges_pos.at(tmp_random_i).erase(tmp_random_j);
        edges_pos.at(tmp_random_j).erase(tmp_random_i);
        edges_neg.at(tmp_random_i).insert(tmp_random_j);
        edges_neg.at(tmp_random_j).insert(tmp_random_i);
      } else {
        // Rcout << "Changed from 0 to -1" << std::endl;
        edges_neg.at(tmp_random_i).insert(tmp_random_j);
        edges_neg.at(tmp_random_j).insert(tmp_random_i);
      }
    }
    // clock.tock("Update");
    
    // clock.tock("Change network");
    
  }
  // clock.stop("clock");
}

void simulate_network_mh( arma::vec coefs,
                          std::unordered_map< int, std::set<int>> &edges_pos,
                          std::unordered_map< int, std::set<int>> &edges_neg,
                          int &n_actors, 
                          int &n_proposals, 
                          int seed, 
                          std::vector<arma::mat> &data_list, 
                          std::vector<int> &type_list,
                          std::vector<arma::vec(*)(const  std::unordered_map< int, std::set<int>> &edges_pos,
                                                const std::unordered_map< int, std::set<int>> &edges_neg,
                                                int &n_actors,
                                                int &tmp_random_i,
                                                int &tmp_random_j, 
                                                arma::mat &data, 
                                                int &type)> &functions, 
                                                arma::vec &global_stats) {
  // Set up objects 
  // Plan is to implement TNT sampling in the future ... 
  // where edges are selected with prob. 0.5
  // std::set<char> nonzero_edges;
  int present_val, proposed_change;
  // double MR;
  arma::mat HR;
  // change_staCt = arma::vec(n_elem);
  arma::mat change_stat;
  set_seed(seed);
  NumericVector random_accept= runif(n_proposals,0,1);
  NumericVector random_propose= runif(n_proposals,0,1);
  arma::vec tmp_vec, random_i(n_proposals), random_j(n_proposals), tmp_stat;
  arma::umat comp_vec;
  random_i = simulate_numbers(0,n_actors-1, n_proposals)+1;
  random_j = simulate_numbers(0,n_actors-1, n_proposals)+1;
  // Where is random_i bigger than random_j? 
  comp_vec = (random_i>random_j);
  arma::uvec ids = find(comp_vec == 1); // Find indices
  tmp_vec = random_i(ids);
  random_i(ids) = random_j(ids);
  random_j(ids) = tmp_vec;
  comp_vec = (random_i==random_j);
  ids = find(comp_vec == 0); // Find indices
  random_i = random_i(ids);
  random_j = random_j(ids);
  int n_proposals_effective = random_j.size();
  // Rcout << n_proposals_effective << std::endl;
  
  // Rcpp::Clock clock;
  int tmp_random_i,tmp_random_j; 
  arma::vec tmp_row;
  // Go through a loop for each proposed change
  for(int a = 0; a <=(n_proposals_effective-1); a = a + 1 ) {
    // clock.tick("Pick present val");
    // Metropolis Hastings or Gibbs Sampling? 
    tmp_random_i = random_i.at(a);
    tmp_random_j = random_j.at(a);
    
    // If the entries are the same (implying loops) we skip the proposal
    // (tmp_random_i==tmp_random_j) ? continue;
    
    
    // Here we pick a random entry
    // present_val = network(tmp_random_i, tmp_random_j);
    if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
      present_val = 1;
    } else if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
      present_val = -1;
    } else {
      present_val = 0;
    }
    // clock.tock("Time to pick");
    // clock.tock("Pick present val");
    // For Gibbs sampling we know that the conditional probability for an 
    // entry to be -1, 0 and 1 is given by the multinomial distribution
    // where each probability is proportional to exp(\theta s(y_ij^(+/-/0)))
    // 1. Step: Calculate for tmp_i[a], tmp_j[a] (i,j) the resulting sufficient statistics under the 
    //          assumption that y_{tmp_i[a], tmp_j[a]} is -1, 0, or 1 
    // 2. Step: Sample once from this multinomial distribution and 
    //          set network(tmp_i[a], tmp_j[a]) to the respective value
    // For now only MH will be implemented 
    // For MH sampling we pick the entry tmp_entry = network(tmp_i[a], tmp_j[a]) 
    //  and depending on its value we have to calculate two change statistics 
    // 1. Step: Check the value of this entry and sample a new proposal 
    // If tmp_entry = 1:
    //    Sample between 0 and -1 with equal probability (stc)
    
    if(present_val == 1){
      if(random_propose(a)>0.5){
        proposed_change = 0;
      } else {
        proposed_change = -1;
      }
    }
    // If present_val = -1:
    //    Sample between 0 and 1 with equal probability (stc)
    if(present_val == -1){
      if(random_propose(a)>0.5){
        proposed_change = 0;
      } else {
        proposed_change = 1;
      }
    }
    // If present_val = 0:
    //    Sample between 1 and -1 with equal probability (stc)
    if(present_val == 0){
      if(random_propose(a)>0.5){
        proposed_change = 1;
      } else {
        proposed_change = -1;
      }
    }
    
    
    // clock.tick("Calculate change stats");
    
    // TNT Sampled would be biased towards 1->0 and -1->0 
    // since most networks are sparse 
    // 2. Step: Calculate the change_stats for i,j for all terms -> we calculate n termsx3 matrix (0 to 1,0 to -1, 1 to -1)
    change_stat = calculate_change_stats(tmp_random_i,
                                         tmp_random_j,
                                         n_actors,
                                         edges_pos,
                                         edges_neg,  
                                         data_list, 
                                         type_list,
                                         functions);
    // This function is only made for changes where present_val and proposed_change are different from each other  
    tmp_row = pick_row(present_val,proposed_change);
    // 3. Step: Calculate the Hastings Ratios by exp(delta(tmp_entry)*theta)
    // Rcout << tmp_row.at(1) << std::endl;
    // Rcout << change_stat << std::endl;
    // 
    // Rcout << change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0) << std::endl;
    // Rcout << change_stat << std::endl;
    // Rcout << tmp_row << std::endl;
    // Rcout << change_stat.col(tmp_row.at(1)-1) << std::endl;
    
    tmp_stat=change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0);
    HR= exp(coefs.t()*tmp_stat);
    // Rcout << "Changed from 1 to 0" << std::endl;
    // Rcout << "Changed from 1 to 0" << std::endl;
    // clock.tock("Calculate change stats");
    
    // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
    // clock.tick("Change val");
    if(random_accept(a)<HR.at(0)){
      // clock.tick("Change network");
      global_stats += tmp_stat;
      // Here we modify the network
      if(proposed_change == 0){
        if(present_val == 1){
          // Rcout << "Changed from 1 to 0" << std::endl;
          edges_pos.at(tmp_random_i).erase(tmp_random_j);
          edges_pos.at(tmp_random_j).erase(tmp_random_i);
        } else {
          // Rcout << "Changed from -1 to 0" << std::endl;
          edges_neg.at(tmp_random_i).erase(tmp_random_j);
          edges_neg.at(tmp_random_j).erase(tmp_random_i);
        }
      }
      if(proposed_change == 1){
        if(present_val == -1){
          // Rcout << "Changed from -1 to 1" << std::endl;
          edges_neg.at(tmp_random_i).erase(tmp_random_j);
          edges_neg.at(tmp_random_j).erase(tmp_random_i);
          edges_pos.at(tmp_random_i).insert(tmp_random_j);
          edges_pos.at(tmp_random_j).insert(tmp_random_i);
        } else {
          // Rcout << "Changed from -1 to 0" << std::endl;
          edges_pos.at(tmp_random_i).insert(tmp_random_j);
          edges_pos.at(tmp_random_j).insert(tmp_random_i);
        }
      }
      if(proposed_change == -1){
        if(present_val == 1){
          // Rcout << "Changed from 1 to -1" << std::endl;
          edges_pos.at(tmp_random_i).erase(tmp_random_j);
          edges_pos.at(tmp_random_j).erase(tmp_random_i);
          edges_neg.at(tmp_random_i).insert(tmp_random_j);
          edges_neg.at(tmp_random_j).insert(tmp_random_i);
        } else {
          // Rcout << "Changed from 0 to -1" << std::endl;
          edges_neg.at(tmp_random_i).insert(tmp_random_j);
          edges_neg.at(tmp_random_j).insert(tmp_random_i);
        }
      }
      // clock.tock("Change network");
      
    }
    // clock.tock("Change val");
  }
  // clock.stop("clock");
  // return(List::create(Named("edges_pos") =edges_pos,
  //                     _["edges_neg"] = edges_neg));
  // network = map_to_mat(edges_pos, edges_neg,n_actors);
  // return network;
}

List simulate_network_alt( arma::vec coefs,
                           std::unordered_map< int, std::set<int>> &edges_pos,
                           std::unordered_map< int, std::set<int>> &edges_neg,
                           int &n_actors, 
                           int &n_proposals, 
                           int seed, 
                           std::vector<arma::mat> &data_list, 
                           std::vector<int> &type_list,
                           std::vector<arma::vec(*)(const  std::unordered_map< int, std::set<int>> &edges_pos,
                                                 const std::unordered_map< int, std::set<int>> &edges_neg,
                                                 int &n_actors,
                                                 int &tmp_random_i,
                                                 int &tmp_random_j, 
                                                 arma::mat &data, 
                                                 int &type)> &functions, 
                                                 arma::vec &global_stats) {
  // Set up objects 
  // Plan is to implement TNT sampling in the future ... 
  // where edges are selected with prob. 0.5
  // std::set<char> nonzero_edges;
  int present_val, proposed_change;
  // double MR;
  arma::mat HR;
  // change_staCt = arma::vec(n_elem);
  arma::mat change_stat;
  set_seed(seed);
  NumericVector random_accept= runif(n_proposals,0,1);
  NumericVector random_propose= runif(n_proposals,0,1);
  arma::vec tmp_vec, random_i(n_proposals), random_j(n_proposals), tmp_stat;
  arma::umat comp_vec;
  random_i = simulate_numbers(0,n_actors-1, n_proposals)+1;
  random_j = simulate_numbers(0,n_actors-1, n_proposals)+1;
  // Where is random_i bigger than random_j? 
  comp_vec = (random_i>random_j);
  arma::uvec ids = find(comp_vec == 1); // Find indices
  tmp_vec = random_i(ids);
  random_i(ids) = random_j(ids);
  random_j(ids) = tmp_vec;
  comp_vec = (random_i==random_j);
  ids = find(comp_vec == 0); // Find indices
  random_i = random_i(ids);
  random_j = random_j(ids);
  int n_proposals_effective = random_j.size();
  // Rcout << n_proposals_effective << std::endl;
  
  // Rcpp::Clock clock;
  int tmp_random_i,tmp_random_j; 
  arma::vec tmp_row;
  std::vector<int> i_val;
  std::vector<int> j_val;
  std::vector<int> change_val;
  int n = 0;
  // Go through a loop for each proposed change
  for(int a = 0; a <=(n_proposals_effective-1); a = a + 1 ) {
    // clock.tick("Pick present val");
    // Metropolis Hastings or Gibbs Sampling? 
    tmp_random_i = random_i.at(a);
    tmp_random_j = random_j.at(a);
    
    // If the entries are the same (implying loops) we skip the proposal
    // (tmp_random_i==tmp_random_j) ? continue;
    
    
    // Here we pick a random entry
    // present_val = network(tmp_random_i, tmp_random_j);
    if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
      present_val = 1;
    } else if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
      present_val = -1;
    } else {
      present_val = 0;
    }
    // clock.tock("Time to pick");
    // clock.tock("Pick present val");
    // For Gibbs sampling we know that the conditional probability for an 
    // entry to be -1, 0 and 1 is given by the multinomial distribution
    // where each probability is proportional to exp(\theta s(y_ij^(+/-/0)))
    // 1. Step: Calculate for tmp_i[a], tmp_j[a] (i,j) the resulting sufficient statistics under the 
    //          assumption that y_{tmp_i[a], tmp_j[a]} is -1, 0, or 1 
    // 2. Step: Sample once from this multinomial distribution and 
    //          set network(tmp_i[a], tmp_j[a]) to the respective value
    // For now only MH will be implemented 
    // For MH sampling we pick the entry tmp_entry = network(tmp_i[a], tmp_j[a]) 
    //  and depending on its value we have to calculate two change statistics 
    // 1. Step: Check the value of this entry and sample a new proposal 
    // If tmp_entry = 1:
    //    Sample between 0 and -1 with equal probability (stc)
    
    if(present_val == 1){
      if(random_propose(a)>0.5){
        proposed_change = 0;
      } else {
        proposed_change = -1;
      }
    }
    // If present_val = -1:
    //    Sample between 0 and 1 with equal probability (stc)
    if(present_val == -1){
      if(random_propose(a)>0.5){
        proposed_change = 0;
      } else {
        proposed_change = 1;
      }
    }
    // If present_val = 0:
    //    Sample between 1 and -1 with equal probability (stc)
    if(present_val == 0){
      if(random_propose(a)>0.5){
        proposed_change = 1;
      } else {
        proposed_change = -1;
      }
    }
    
    
    // clock.tick("Calculate change stats");
    
    // TNT Sampled would be biased towards 1->0 and -1->0 
    // since most networks are sparse 
    // 2. Step: Calculate the change_stats for i,j for all terms -> we calculate ntermsx3 matrix (0 to 1,0 to -1, 1 to -1)
    change_stat = calculate_change_stats(tmp_random_i,
                                         tmp_random_j,
                                         n_actors,
                                         edges_pos,
                                         edges_neg,  
                                         data_list, 
                                         type_list,
                                         functions);
    // This function is only made for changes where present_val and proposed_change are different from each other  
    tmp_row = pick_row(present_val,proposed_change);
    // 3. Step: Calculate the Hastings Ratios by exp(delta(tmp_entry)*theta)
    // Rcout << tmp_row.at(1) << std::endl;
    // Rcout << change_stat << std::endl;
    // 
    // Rcout << change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0) << std::endl;
    // Rcout << change_stat << std::endl;
    // Rcout << tmp_row << std::endl;
    // Rcout << change_stat.col(tmp_row.at(1)-1) << std::endl;
    
    tmp_stat=change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0);
    HR= exp(coefs.t()*tmp_stat);
    // Rcout << "Changed from 1 to 0" << std::endl;
    // Rcout << "Changed from 1 to 0" << std::endl;
    // clock.tock("Calculate change stats");
    
    // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
    // clock.tick("Change val");
    if(random_accept(a)<HR.at(0)){
      i_val.push_back(tmp_random_i);
      j_val.push_back(tmp_random_j);
      change_val.push_back(proposed_change);
      n +=1;
      // clock.tick("Change network");
      global_stats += tmp_stat;
      // Here we modify the network
      if(proposed_change == 0){
        if(present_val == 1){
          // Rcout << "Changed from 1 to 0" << std::endl;
          edges_pos.at(tmp_random_i).erase(tmp_random_j);
          edges_pos.at(tmp_random_j).erase(tmp_random_i);
        } else {
          // Rcout << "Changed from -1 to 0" << std::endl;
          edges_neg.at(tmp_random_i).erase(tmp_random_j);
          edges_neg.at(tmp_random_j).erase(tmp_random_i);
        }
      }
      if(proposed_change == 1){
        if(present_val == -1){
          // Rcout << "Changed from -1 to 1" << std::endl;
          edges_neg.at(tmp_random_i).erase(tmp_random_j);
          edges_neg.at(tmp_random_j).erase(tmp_random_i);
          edges_pos.at(tmp_random_i).insert(tmp_random_j);
          edges_pos.at(tmp_random_j).insert(tmp_random_i);
        } else {
          // Rcout << "Changed from -1 to 0" << std::endl;
          edges_pos.at(tmp_random_i).insert(tmp_random_j);
          edges_pos.at(tmp_random_j).insert(tmp_random_i);
        }
      }
      if(proposed_change == -1){
        if(present_val == 1){
          // Rcout << "Changed from 1 to -1" << std::endl;
          edges_pos.at(tmp_random_i).erase(tmp_random_j);
          edges_pos.at(tmp_random_j).erase(tmp_random_i);
          edges_neg.at(tmp_random_i).insert(tmp_random_j);
          edges_neg.at(tmp_random_j).insert(tmp_random_i);
        } else {
          // Rcout << "Changed from 0 to -1" << std::endl;
          edges_neg.at(tmp_random_i).insert(tmp_random_j);
          edges_neg.at(tmp_random_j).insert(tmp_random_i);
        }
      }
      // clock.tock("Change network");
      
    }
    // clock.tock("Change val");
  }
  // clock.stop("clock");
  // return(List::create(Named("edges_pos") =edges_pos,
  //                     _["edges_neg"] = edges_neg));
  // network = map_to_mat(edges_pos, edges_neg,n_actors);
  // return network;
  return(List::create(Named("i_vals") = i_val , _["j_vals"] = j_val, _["change_vals"] = change_val));
  
}

// // Count global statistics 
// arma::vec count_statistics_new(std::unordered_map< int, std::set<int>> &edges_pos,
//                                std::unordered_map< int, std::set<int>> &edges_neg,
//                                int n_actors){
//   arma::vec res(2);
//   edges_pos.bucket_size(1);
//   for(int a = 1; a <=(edges_pos.size()); a = a + 1 ) {
//     res.at(0) += edges_pos.at(a).size();
//     res.at(1) += edges_neg.at(a).size();
//   }
//   return(res);
// }


// [[Rcpp::export]]
List simulate_networks_new( arma::vec& coefs, 
                            std::vector<std::string>& terms,
                            int& n_actors,
                            std::vector<arma::mat>& data_list, 
                            std::vector<int>& type_list,
                            bool mh,
                            int n_proposals = 100,
                            int n_proposals_burn_in = 100, 
                            int seed = 123, 
                            int number_networks = 1, 
                            bool only_stats = false){
  std::unordered_map< int, std::set<int> > edges_pos;
  for (int i = 1; i <= n_actors; i++){
    edges_pos[i] = std::set<int>();
  }
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  List res; 
  
  arma::mat stats(number_networks,terms.size());
  stats.fill(0);
  // Generate change statistic function from the terms 
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  
  functions = change_statistics_generate(terms);
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  // Start for a burn in period with the normal number of proposals 
  // Intialize global statistics and then adapt them peu a peu 
  arma::vec global_stats(functions.size());
  global_stats.fill(0); 
  global_stats = global_stats + at_zero;
  if(mh){
    simulate_network_mh(coefs,edges_pos, edges_neg, n_actors, 
                        n_proposals_burn_in, seed,
                        data_list, type_list,functions, 
                        global_stats);
  } else { 
    simulate_network_gs(coefs,edges_pos, edges_neg, n_actors, 
                        n_proposals_burn_in, seed,
                        data_list, type_list,functions, 
                        global_stats);
  }
  for(int i = 0; i <(number_networks);i ++) {
    // Simulate the network
    if(mh){
      simulate_network_mh(coefs,edges_pos, edges_neg, n_actors, 
                          n_proposals, seed +i,
                          data_list, type_list,functions, 
                          global_stats);
    } else {
      simulate_network_gs(coefs,edges_pos, edges_neg, n_actors, 
                          n_proposals, seed + i,
                          data_list, type_list,functions, 
                          global_stats);
    }
    if(only_stats){
      // Count global statistics 
      stats.row(i) = global_stats.as_row();
    } else{
      // Save the networks
      res.push_back(List::create(Named("edges_pos") = edges_pos , _["edges_neg"] = edges_neg));
      // Count global statistics 
      stats.row(i) = global_stats.as_row();
    }

  }
  if(only_stats){
    return(List::create(_["stats"] = stats));
  } else {
    return(List::create(Named("networks") = res , _["stats"] = stats));
  }
  // clock.stop("clock");
}
// [[Rcpp::export]]
arma::mat simulate_networks_stats(arma::vec& coefs, 
                            std::vector<std::string>& terms,
                            int& n_actors,
                            std::vector<arma::mat>& data_list, 
                            std::vector<int>& type_list,
                            bool mh,
                            int n_proposals = 100,
                            int n_proposals_burn_in = 100, 
                            int seed = 123, 
                            int number_networks = 1){
  std::unordered_map< int, std::set<int> > edges_pos;
  for (int i = 1; i <= n_actors; i++){
    edges_pos[i] = std::set<int>();
  }
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;

  arma::mat stats(number_networks,terms.size());
  stats.fill(0);
  // Generate change statistic function from the terms 
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  
  functions = change_statistics_generate(terms);
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  // Start for a burn in period with the normal number of proposals 
  // Intialize global statistics and then adapt them peu a peu 
  arma::vec global_stats(functions.size());
  global_stats.fill(0); 
  global_stats = global_stats + at_zero;
  if(mh){
    simulate_network_mh(coefs,edges_pos, edges_neg, n_actors, 
                        n_proposals_burn_in, seed,
                        data_list, type_list,functions, 
                        global_stats);
  } else {  
    simulate_network_gs(coefs,edges_pos, edges_neg, n_actors, 
                        n_proposals_burn_in, seed,
                        data_list, type_list,functions, 
                        global_stats);
  } 
  for(int i = 0; i <(number_networks);i ++) {
    // Simulate the network
    if(mh){
      simulate_network_mh(coefs,edges_pos, edges_neg, n_actors, 
                          n_proposals, seed +i,
                          data_list, type_list,functions, 
                          global_stats);
    } else { 
      simulate_network_gs(coefs,edges_pos, edges_neg, n_actors, 
                          n_proposals, seed + i,
                          data_list, type_list,functions, 
                          global_stats);
    } 
    stats.row(i) = global_stats.as_row();
  }
  return(stats);
  // clock.stop("clock");
} 


// [[Rcpp::export]]
List simulate_trying( arma::vec& coefs, 
                      std::vector<std::string>& terms,
                      int& n_actors,
                      std::vector<arma::mat>& data_list, 
                      std::vector<int>& type_list,
                      int n_proposals = 100,
                      int n_proposals_burn_in = 100, 
                      int seed = 123, 
                      int number_networks = 1, 
                      bool only_stats = false){
  std::unordered_map< int, std::set<int> > edges_pos;
  for (int i = 1; i <= n_actors; i++){
    edges_pos[i] = std::set<int>();
  }
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  
  arma::mat stats(number_networks,terms.size());
  stats.fill(0);
  // Generate change statistic function from the terms 
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  
  functions = change_statistics_generate(terms);
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  // Start for a burn in period with the normal number of proposals 
  // Intialize global statistics and then adapt them peu a peu 
  arma::vec global_stats(functions.size());
  global_stats.fill(0); 
  global_stats = global_stats + at_zero;
  List res;
  res.push_back(simulate_network_alt(coefs,edges_pos, edges_neg, n_actors, 
                                     n_proposals_burn_in, seed,
                                     data_list, type_list,functions, 
                                     global_stats));
  ;
  for(int i = 0; i <(number_networks);i ++) {
    // Simulate the network
    res.push_back(simulate_network_alt(coefs,edges_pos, edges_neg, n_actors, 
                                       n_proposals, seed +i,
                                       data_list, type_list,functions, 
                                       global_stats));
    stats.row(i) = global_stats.as_row();
  }
  return(List::create(Named("res") = res , _["stats"] = stats));
}




arma::mat simulate_networks_intern(arma::vec& coefs, 
                                   std::unordered_map< int, std::set<int>>& edges_pos, 
                                   std::unordered_map< int, std::set<int>>& edges_neg,
                                   arma::vec global_stats,
                                   bool mh,
                                   std::vector<arma::vec(*)(const  std::unordered_map< int, std::set<int>> &edges_pos,
                                                         const std::unordered_map< int, std::set<int>> &edges_neg,
                                                         int &n_actors,
                                                         int &present_val,
                                                         int &proposed_change, 
                                                         arma::mat &data, 
                                                         int &type)> functions,
                                                         int& n_actors,
                                                         std::vector<arma::mat>& data_list, 
                                                         std::vector<int>& type_list,
                                                         int n_proposals = 100,
                                                         int n_proposals_burn_in = 100, 
                                                         int seed = 123, 
                                                         int number_networks = 1){
  
  arma::mat stats(number_networks,functions.size());
  stats.fill(0);
  // Start for a burn in period with the normal number of proposals 
  // Intialize global statistics and then adapt them peu a peu 
  // arma::vec global_stats(functions.size());
  // global_stats.fill(0);
  
  
  if(mh){
    simulate_network_mh(coefs,edges_pos, edges_neg, n_actors, 
                        n_proposals_burn_in, seed,
                        data_list, type_list,functions, 
                        global_stats);
  } else {
    simulate_network_gs(coefs,edges_pos, edges_neg, n_actors, 
                        n_proposals_burn_in, seed,
                        data_list, type_list,functions, 
                        global_stats);
  }
  
  for(int i = 0; i <(number_networks);i ++) {
    // Simulate the network
    if(mh){
      simulate_network_mh(coefs,edges_pos, edges_neg, n_actors, 
                          n_proposals, seed +i,
                          data_list, type_list,functions, 
                          global_stats);
    } else {
      // Rcout << "Changed from 1 to 0" << std::endl;
      simulate_network_gs(coefs,edges_pos, edges_neg, n_actors, 
                          n_proposals, seed+i,
                          data_list, type_list,functions, 
                          global_stats);
    }
    // Count global statistics 
    // Rcout << "Changed from 1 to 0" << std::endl;
    stats.row(i) = global_stats.as_row();
  }
  return(stats);
  // clock.stop("clock");
}


// Count global statistics
arma::vec count_global_statistic( std::unordered_map< int, std::set<int>> &edges_pos,
                                  std::unordered_map< int, std::set<int>> &edges_neg,
                                  int n_actors, 
                                  std::vector<arma::mat> &data_list, 
                                  std::vector<int> &type_list,
                                  std::vector<arma::vec(*)(const  std::unordered_map< int, std::set<int>> &edges_pos,
                                                        const  std::unordered_map< int, std::set<int>> &edges_neg,
                                                        int &n_actors,
                                                        int &tmp_random_i,
                                                        int &tmp_random_j, 
                                                        arma::mat &data, 
                                                        int &type)> functions) {
  
  std::unordered_map< int, std::set<int>> edges_pos_alt;
  for (int i = 1; i <= n_actors; i++){
    edges_pos_alt[i] = std::set<int>();
  }
  std::unordered_map< int, std::set<int>> edges_neg_alt = edges_pos_alt;
  
  arma::vec res(functions.size());
  arma::mat change_stat;
  arma::vec tmp_row;
  std::set<int> tmp_js;
  for (int i = 1; i <= n_actors; i++){
    // Rcout << "Start with" << std::endl;
    tmp_js = edges_pos.at(i);
    // std::vector<int> tmp_js_vector(tmp_js.begin(), tmp_js.end());
    // std::vector<size_t> y(tmp_js_vector.size());
    // 
    // std::iota(y.begin(), y.end(), 0);
    // std::copy_if(y.begin(), y.end(),
    //              std::ostream_iterator<size_t>(std::cout, " "),
    //              [&](size_t j) { return tmp_js_vector[j] > i; });
    // 
    // 
    // Rcout << "tmp_js_vector" << std::endl;
    // for (int x: tmp_js_vector){
    //   Rcout << x << std::endl;
    // }
    // 
    // Rcout << "tmp_js" << std::endl;
    // for (int x: tmp_js) {
    //   Rcout << x << std::endl;
    // }
    // 
    if(tmp_js.size()>0){
      std::set<int>::iterator it = tmp_js.begin();
      
      while (it != tmp_js.end())
      {
        if(*it>i){
          change_stat = calculate_change_stats(i, 
                                               *it,
                                               n_actors,
                                               edges_pos_alt,
                                               edges_neg_alt,
                                               data_list, 
                                               type_list,
                                               functions);
          edges_pos_alt.at(i).insert(*it);
          edges_pos_alt.at(*it).insert(i);
          tmp_row = pick_row(0,1);
          // Rcout << "Adding Positive" << std::endl; 
          // Rcout << i << std::endl;
          // Rcout << *it << std::endl;
          // Rcout << change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0) << std::endl;
          // 
          res +=change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0); 
        }
        it++;
      }
    }
    
    tmp_js = edges_neg.at(i);
    if(tmp_js.size()>0){
      std::set<int>::iterator it = tmp_js.begin();
      
      while (it != tmp_js.end())
      {
        if(*it>i){
          
          change_stat = calculate_change_stats(i, 
                                               *it,
                                               n_actors,
                                               edges_pos_alt,
                                               edges_neg_alt,
                                               data_list, 
                                               type_list,
                                               functions);
          
          tmp_row = pick_row(0,-1);
          edges_neg_alt.at(i).insert(*it);
          edges_neg_alt.at(*it).insert(i);
          // Rcout << "Adding Negative" << std::endl; 
          // Rcout << i << std::endl;
          // Rcout << *it << std::endl;
          // Rcout << change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0) << std::endl;
          
          res +=change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0);
        }
        it++;
      }
    }
    
    
  }
  return(res);
}
// [[Rcpp::export]]
arma::vec count_global(arma::mat network ,std::vector<std::string> terms, int n_actors,
                       std::vector<arma::mat> &data_list, 
                       std::vector<int> &type_list) { 
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  List coefs, tmp_samples;
  arma::mat tmp_mean,tmp_stats, tmp_var;
  arma::vec update_coef,tmp_coef;
  NumericVector beg_coef;
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  arma::vec at_zero; 
  at_zero = eval_at_empty_network(terms, n_actors);
  
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg, 
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  return(global_stats);
  
}

// [[Rcpp::export]]
List trying_again2(arma::mat network ,std::vector<std::string> terms, int n_actors,
                   std::vector<arma::mat> &data_list, 
                   std::vector<int> &type_list, 
                   std::vector<int> from, 
                   std::vector<int> to, 
                   std::vector<int> change) { 
  
  
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  List coefs, tmp_samples;
  arma::mat tmp_mean,tmp_stats, tmp_var;
  arma::vec update_coef,tmp_coef;
  NumericVector beg_coef;
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  arma::vec at_zero; 
  at_zero = eval_at_empty_network(terms, n_actors);
  
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg, 
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  // global_stats = global_stats + at_zero;
  Rcout << "Starding stats are: "<< std::endl;
  Rcout << global_stats<< std::endl;
  arma::mat change_stat;
  arma::vec tmp_row;
  for (unsigned int j = 0; j < change.size(); j++){
    change_stat = calculate_change_stats(from.at(j), 
                                         to.at(j),
                                         n_actors,
                                         edges_pos,
                                         edges_neg,
                                         data_list, 
                                         type_list,
                                         functions);
    tmp_row = pick_row(network.at(from.at(j) -1, to.at(j)-1),change.at(j));
    Rcout << "Change at " + std::to_string(from.at(j)) + " and " + std::to_string(to.at(j))<< std::endl;
    Rcout << "From " + std::to_string(network.at(from.at(j) -1, to.at(j)-1)) + " to " + std::to_string(change.at(j))<< std::endl;
    
    if(network.at(from.at(j) -1, to.at(j)-1) != change.at(j)){
      Rcout << "Change is "<< std::endl;
      Rcout << change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0)<< std::endl;
      global_stats +=change_stat.col(tmp_row.at(1)-1)*tmp_row.at(0); 
      if(change.at(j) == -1){
        if(network.at(from.at(j) -1, to.at(j)-1) == 0){
          edges_neg.at(from.at(j)).insert(to.at(j));
          edges_neg.at(to.at(j)).insert(from.at(j));
        } else {
          edges_neg.at(from.at(j)).insert(to.at(j));
          edges_neg.at(to.at(j)).insert(from.at(j));
          edges_pos.at(from.at(j)).erase(to.at(j));
          edges_pos.at(to.at(j)).erase(from.at(j));
        }
      }
      if(change.at(j) == 1){
        if(network.at(from.at(j) -1, to.at(j)-1) == 0){
          edges_pos.at(from.at(j)).insert(to.at(j));
          edges_pos.at(to.at(j)).insert(from.at(j));
        } else {
          edges_pos.at(from.at(j)).insert(to.at(j));
          edges_pos.at(to.at(j)).insert(from.at(j));
          edges_neg.at(from.at(j)).erase(to.at(j));
          edges_neg.at(to.at(j)).erase(from.at(j));
        } 
      } 
      if(change.at(j) == 0){
        if(network.at(from.at(j) -1, to.at(j)-1) == 1){
          edges_pos.at(from.at(j)).erase(to.at(j));
          edges_pos.at(to.at(j)).erase(from.at(j));
        } else {
          edges_neg.at(from.at(j)).erase(to.at(j));
          edges_neg.at(to.at(j)).erase(from.at(j));
        }
      }  
      network.at(from.at(j) -1, to.at(j)-1) = change.at(j);
      network.at(to.at(j) -1, from.at(j)-1) = change.at(j);
    } else{
      Rcout << "No Change"<< std::endl;
    }
    
    Rcout << "Leading to "<< std::endl;
    Rcout << global_stats<< std::endl;
  }
  return(List::create(Named("mat") =network,
                      _["stats"] = global_stats));
  
}

// Function to check if the simulation starting on the observed network works 
// [[Rcpp::export]]
List simulation_mat(arma::mat network ,
                  std::vector<std::string> terms, 
                  int n_actors,
                  int n_proposals,
                  int seed,
                  int number_networks, 
                  std::vector<arma::mat> &data_list, 
                  std::vector<int> &type_list,
                  arma::vec &coef, 
                  bool mh) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  List res;
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg,
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  res.push_back(global_stats);
  arma::mat tmp_stats;
  
  tmp_stats = simulate_networks_intern(coef,edges_pos, edges_neg,global_stats, mh,functions, n_actors,
                                       data_list, type_list,
                                       n_proposals,n_proposals, 
                                       seed,number_networks);
  res.push_back(tmp_stats);
  return(res);
}  
// [[Rcpp::export]]
arma::vec count_global_terms(arma::mat network ,
                             std::vector<std::string> terms, 
                             int n_actors,
                             std::vector<arma::mat> &data_list, 
                             std::vector<int> &type_list) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const  std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg,
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  return(global_stats);
}

// Here the mle estimation is done  according to Hunter    
// [[Rcpp::export]]     
std::vector<arma::vec> mle_estimation(arma::mat network ,
                                      std::vector<std::string> terms, 
                                      int n_actors,
                                      int n_proposals,
                                      int n_proposals_burn_in,
                                      int seed,
                                      int number_networks, 
                                      std::vector<arma::mat> &data_list, 
                                      std::vector<int> &type_list,
                                      arma::vec &beg_coef,
                                      bool mh,
                                      int max_it = 10, 
                                      double tol = 0.001, 
                                      bool start_with_empty_net = true) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean,tmp_stats, tmp_var;
  arma::vec update_coef,tmp_coef; 
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  

  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg,
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  // Set \theta_0
  // Rcout << global_stats << std::endl;
  coefs.push_back(beg_coef);
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl; 
    tmp_coef = coefs.at(i);
    if(start_with_empty_net){
      tmp_stats = simulate_networks_stats(tmp_coef,terms, n_actors,
                                          data_list, type_list,mh,
                                          n_proposals,n_proposals_burn_in, 
                                          seed +i,number_networks);
    } else { 
      tmp_stats = simulate_networks_intern(tmp_coef,
                                           edges_pos, 
                                           edges_neg,
                                           global_stats,
                                           mh, 
                                           functions, 
                                           n_actors,
                                           data_list, 
                                           type_list,
                                           n_proposals,
                                           n_proposals_burn_in,
                                           seed +i,
                                           number_networks);
    }
    
    // Rcout << tmp_stats << std::endl;
    // Rcout << tmp_coef << std::endl;
    // tmp_stats =as<arma::mat> (tmp_samples["stats"]);
    // Rcout << tmp_stats << std::endl;
    tmp_var = arma::cov(tmp_stats);
    // Rcout << tmp_var << std::endl;
    
    // Rcout << tmp_var << std::endl;
    tmp_mean = arma::mean(tmp_stats,0);
    update_coef = tmp_coef + arma::solve(tmp_var,global_stats-tmp_mean.t());
    
    // update_coef = tmp_coef + arma::inv(tmp_var)*(global_stats-tmp_mean.t());
    coefs.push_back(update_coef);
    // Rcout << tmp_mean << std::endl;
    // Rcout << update_coef << std::endl;
    
    bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    if(is_converged){
      break;
    }
    // coefs.push_back(arma::var(tmp_stats));
  }
  return(coefs);
}

// [[Rcpp::export]]     
std::vector<arma::vec> mle_estimation_stepping(arma::mat network ,
                                               std::vector<std::string> terms, 
                                               int n_actors,
                                               int n_proposals,
                                               int n_proposals_burn_in,
                                               int seed,
                                               int number_networks, 
                                               std::vector<arma::mat> &data_list, 
                                               std::vector<int> &type_list,
                                               arma::vec &beg_coef,
                                               int steps,
                                               bool mh,
                                               int max_it = 10, 
                                               double tol = 0.001, 
                                               bool start_with_empty_net = true) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean,tmp_stats, tmp_var;
  arma::vec update_coef,tmp_coef; 
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  

  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg,
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  // Set \theta_0
  // Rcout << global_stats << std::endl;
  coefs.push_back(beg_coef);
  // Generate the sequence of possible gamma values
  arma::vec poss_gamma = arma::linspace(0, 1, steps);
  // Rcout << poss_gamma << std::endl;
  NumericVector n;
  NumericVector one = {1};
  NumericVector old_gamma = {0};;
  // This function is taken from ergm, which in turn uses lpSolveAPI
  Function f("find_min_in_ch");   
  
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl; 
    tmp_coef = coefs.at(i);
    if(start_with_empty_net){
      tmp_stats = simulate_networks_stats(tmp_coef,terms, n_actors,
                                          data_list, type_list,mh,
                                          n_proposals,n_proposals_burn_in, 
                                          seed +i,number_networks);
    } else { 
      tmp_stats = simulate_networks_intern(tmp_coef,
                                           edges_pos, 
                                           edges_neg,
                                           global_stats,
                                           mh, 
                                           functions, 
                                           n_actors,
                                           data_list, 
                                           type_list,
                                           n_proposals,
                                           n_proposals_burn_in,
                                           seed +i,
                                           number_networks);
    }
    
    // Rcout << tmp_stats << std::endl;
    // Rcout << tmp_coef << std::endl;
    // tmp_stats =as<arma::mat> (tmp_samples["stats"]);
    // Rcout << tmp_stats << std::endl;
    tmp_var = arma::cov(tmp_stats);
    // Rcout << tmp_var << std::endl;
    // Rcout << tmp_var << std::endl;
    tmp_mean = arma::mean(tmp_stats,0);
    n = f(poss_gamma, global_stats, tmp_stats, tmp_mean);
    Rcout << "Trying gamma="+ std::to_string(poss_gamma.at(n.at(0))) << std::endl;
    
    tmp_global_stats = poss_gamma.at(n.at(0))*global_stats + (1-poss_gamma.at(n.at(0)))*tmp_mean.t();
    update_coef = tmp_coef + arma::solve(tmp_var,tmp_global_stats-tmp_mean.t());
    coefs.push_back(update_coef);
    
    if((old_gamma.at(0) ==  n.at(0)) & ( n.at(0) == (steps-1))){
      break;
    }
    old_gamma = n;
    // Rcout << tmp_mean << std::endl;
    // Rcout << tmp_global_stats << std::endl;
    // Rcout << steps.at(i) << std::endl;
    
    // bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    // if(is_converged){
    //   break;
    // }
    
  }
  // Rcout << "Successfully found the starting value" << std::endl; 
  Rcout << "Start MLE" << std::endl; 
  return(coefs);
}

// [[Rcpp::export]]     
std::vector<arma::vec> t_mle_estimation_stepping(std::vector<arma::mat> networks,
                                               std::vector<std::string> terms, 
                                               int n_actors,
                                               int n_proposals,
                                               int n_proposals_burn_in,
                                               int seed,
                                               int number_networks, 
                                               std::vector<std::vector<arma::mat>> &data_lists, 
                                               std::vector<int> &type_list,
                                               arma::vec &beg_coef,
                                               int steps,
                                               bool mh,
                                               int max_it = 10, 
                                               double tol = 0.001, 
                                               bool start_with_empty_net = true) {
  
  // Set up objects
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_pos;
  std::unordered_map< int, std::set<int>> pos_edges;
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_neg = all_edges_pos;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    all_edges_pos.push_back(mat_to_map(networks.at(i),1, n_actors));
    all_edges_neg.push_back(mat_to_map(networks.at(i),-1, n_actors));
  }
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean(1,terms.size()), tmp_var(terms.size(),terms.size());
  arma::vec update_coef,tmp_coef; 
  std::vector<arma::mat> all_stats;
  arma::mat trying,tmp_stats(number_networks,terms.size()), sum_stats(number_networks,terms.size());
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::mat global_stats_per_network(data_lists.size(),type_list.size());
  arma::vec sum_global_stats(type_list.size()), tmp(type_list.size());
  // Rcout << all_edges_pos.size() << std::endl;
  // Rcout << all_edges_neg.size() << std::endl;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    tmp = count_global_statistic(all_edges_pos.at(i),
                                 all_edges_neg.at(i),
                                 n_actors, 
                                 data_lists.at(i), 
                                 type_list,
                                 functions) ;
    // Rcout << tmp << std::endl;
    // Rcout << at_zero << std::endl;
    // Rcout << tmp + at_zero << std::endl;
    
    tmp += at_zero;
    global_stats_per_network.row(i) = tmp.as_row();
    sum_global_stats += tmp;
  }
  // sum_global_stats = sum_global_stats;
  // global_stats_per_network = global_stats_per_network + at_zero;
  // Rcout << global_stats_per_network << std::endl;
  
  // Rcout << global_stats << std::endl;
  coefs.push_back(beg_coef);
  // Generate the sequence of possible gamma values
  arma::vec poss_gamma = arma::linspace(0, 1, steps);
  // Rcout << poss_gamma << std::endl;
  NumericVector n;
  NumericVector old_gamma = {0};;
  // This function is taken from ergm, which, in turn, uses lpSolveAPI
  Function f("find_min_in_ch");  
  // Function parallel_sample("parallel_sample");
  
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl;
    tmp_coef = coefs.at(i);
    // parallel_sample(tmp_coef,terms, n_actors,
    //                 data_lists, type_list,mh,
    //                 n_proposals,n_proposals_burn_in,
    //                 seed +i,number_networks,cl);
    // for each network simulate the statistics and estimate their mean and covariance 
    for(unsigned int j = 0;j< data_lists.size();j +=1){
      if(start_with_empty_net){
        // a) Sample
        tmp_stats = simulate_networks_stats(tmp_coef,terms, n_actors,
                                                    data_lists.at(j), type_list,mh,
                                                    n_proposals,n_proposals_burn_in,
                                                    seed +i,number_networks);
        all_stats.push_back(tmp_stats);
        sum_stats += tmp_stats;
        // Rcout << sum_stats << std::endl;
        // Rcout << arma::mean(tmp_stats,0) << std::endl;
        // b) Estimate mean
        tmp_mean += arma::mean(tmp_stats,0);
        // c) Estimate variance
        // Rcout << arma::cov(tmp_stats); << std::endl;
        tmp_var += arma::cov(tmp_stats);
      } else {
        tmp_stats=simulate_networks_intern(tmp_coef,
                                           all_edges_pos.at(j),
                                           all_edges_neg.at(j),
                                           global_stats_per_network.row(j).t(),
                                           mh,
                                           functions,
                                           n_actors,
                                           data_lists.at(j),
                                           type_list,
                                           n_proposals,
                                           n_proposals_burn_in,
                                           seed +i,
                                           number_networks);
        // Rcout << trying.size() << std::endl;
        // Rcout << tmp_stats.size() << std::endl;
        // tmp_stats = trying;
        all_stats.push_back(tmp_stats);
        // Rcout << "Done" << std::endl;
        sum_stats += tmp_stats;
        // Rcout << "Done" << std::endl;
        // b) Estimate mean
        tmp_mean += arma::mean(tmp_stats,0);
        // Rcout << "Done" << std::endl;
        // c) Estimate variance
        tmp_var += arma::cov(tmp_stats);
      }
    }
    // Rcout << tmp_mean << std::endl;
    // Rcout << sum_global_stats << std::endl;
    // tmp_var= arma::cov(sum_stats);
    n = f(poss_gamma, sum_global_stats, sum_stats, tmp_mean);
    // Rcout << n.at(0) << std::endl;
    // Rcout <<steps-1 << std::endl;
    Rcout << "Trying gamma="+ std::to_string(poss_gamma.at(n.at(0))) << std::endl;
    // n.at(0) = 1;
    tmp_global_stats = poss_gamma.at(n.at(0))*sum_global_stats + (1-poss_gamma.at(n.at(0)))*tmp_mean.t();
    // tmp_global_stats = sum_global_stats;
    // tmp_global_stats = sum_global_stats/data_lists.size();
    // tmp_var = tmp_var/data_lists.size();
    // tmp_mean = tmp_mean/data_lists.size();
    
    
    update_coef = tmp_coef + arma::solve(tmp_var,tmp_global_stats-tmp_mean.t());
    coefs.push_back(update_coef);

    // Set all values back to zero 
    tmp_var.zeros();
    tmp_mean.zeros();
    sum_stats.zeros();
    
    if((old_gamma.at(0) ==  n.at(0)) & ( n.at(0) == (steps-1))){
      break;
    }
    old_gamma = n;
    // Rcout << tmp_mean << std::endl;
    // Rcout << tmp_global_stats << std::endl;
    // Rcout << steps.at(i) << std::endl;

    // bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    // if(is_converged){
    //   break;
    // }

  }
  Rcout << "Successfully found the starting value" << std::endl;
  Rcout << "Start MLE" << std::endl; 
  return(coefs);
}

// [[Rcpp::export]]     
std::vector<arma::vec> t_mle_estimation_stepping_p(std::vector<arma::mat> networks,
                                                 std::vector<std::string> terms, 
                                                 int n_actors,
                                                 int n_proposals,
                                                 int n_proposals_burn_in,
                                                 int seed,
                                                 int number_networks, 
                                                 std::vector<std::vector<arma::mat>> &data_lists, 
                                                 std::vector<int> &type_list,
                                                 arma::vec &beg_coef,
                                                 int steps,
                                                 bool mh,
                                                 List cluster,
                                                 List data_lists_par, 
                                                 List networks_par, 
                                                 int max_it = 10, 
                                                 double tol = 0.001, 
                                                 bool start_with_empty_net = true) {
  
  // Set up objects
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_pos;
  std::unordered_map< int, std::set<int>> pos_edges;
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_neg = all_edges_pos;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    all_edges_pos.push_back(mat_to_map(networks.at(i),1, n_actors));
    all_edges_neg.push_back(mat_to_map(networks.at(i),-1, n_actors));
  }
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean(1,terms.size()), tmp_var(terms.size(),terms.size());
  arma::vec update_coef,tmp_coef; 
  std::vector<arma::mat> all_stats;
  arma::mat trying,tmp_stats(number_networks,terms.size()), sum_stats(number_networks,terms.size());
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::mat global_stats_per_network(data_lists.size(),type_list.size());
  arma::vec sum_global_stats(type_list.size()), tmp(type_list.size());
  // Rcout << all_edges_pos.size() << std::endl;
  // Rcout << all_edges_neg.size() << std::endl;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    tmp = count_global_statistic(all_edges_pos.at(i),
                                 all_edges_neg.at(i),
                                 n_actors, 
                                 data_lists.at(i), 
                                 type_list,
                                 functions) ;
    // Rcout << tmp << std::endl;
    // Rcout << at_zero << std::endl;
    // Rcout << tmp + at_zero << std::endl;
    
    tmp += at_zero;
    global_stats_per_network.row(i) = tmp.as_row();
    sum_global_stats += tmp;
  }
  // sum_global_stats = sum_global_stats;
  // global_stats_per_network = global_stats_per_network + at_zero;
  // Rcout << global_stats_per_network << std::endl;
  
  // Rcout << global_stats << std::endl;
  coefs.push_back(beg_coef);
  // Generate the sequence of possible gamma values
  arma::vec poss_gamma = arma::linspace(0, 1, steps);
  // Rcout << poss_gamma << std::endl;
  NumericVector n;
  NumericVector old_gamma = {0};;
  // This function is taken from ergm, which, in turn, uses lpSolveAPI
  Function f("find_min_in_ch");  
  Function parallel_sample("parallel_sample");
  List res;
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl;
    tmp_coef = coefs.at(i);
   
   // for each network simulate the statistics and estimate their mean and covariance 
    res = parallel_sample(tmp_coef,terms, n_actors,
                          networks_par,data_lists_par,type_list,mh,
                    n_proposals,n_proposals_burn_in,
                    seed +i,number_networks,cluster,start_with_empty_net);
    // Rcout << n.at(0) << std::endl;
    NumericVector tmp_mean_1 = res.at(0);
     
    NumericMatrix tmp_var_1 = res.at(1);
    NumericMatrix sum_global_stats_1 = res.at(2);
    tmp_mean = arma::mat(tmp_mean_1.begin(),  1, tmp_mean_1.size(),false);
    tmp_var = arma::mat(tmp_var_1.begin(), tmp_var_1.nrow(), tmp_var_1.ncol(), false);
    sum_stats = arma::mat(sum_global_stats_1.begin(), sum_global_stats_1.nrow(), sum_global_stats_1.ncol(), false);
    
    n = f(poss_gamma, sum_global_stats, sum_stats, tmp_mean);
    // Rcout << n.at(0) << std::endl;
    // Rcout <<steps-1 << std::endl;
    Rcout << "Trying gamma="+ std::to_string(poss_gamma.at(n.at(0))) << std::endl;
    // n.at(0) = 1;
    tmp_global_stats = poss_gamma.at(n.at(0))*sum_global_stats + (1-poss_gamma.at(n.at(0)))*tmp_mean.t();
    // tmp_global_stats = sum_global_stats;
    // tmp_global_stats = sum_global_stats/data_lists.size();
    // tmp_var = tmp_var/data_lists.size();
    // tmp_mean = tmp_mean/data_lists.size();
    
    
    update_coef = tmp_coef + arma::solve(tmp_var,tmp_global_stats-tmp_mean.t());
    coefs.push_back(update_coef);
    
    // Set all values back to zero 
    tmp_var.zeros();
    tmp_mean.zeros();
    sum_stats.zeros();
    
    if((old_gamma.at(0) ==  n.at(0)) & ( n.at(0) == (steps-1))){
      break;
    }
    old_gamma = n;
    // Rcout << tmp_mean << std::endl;
    // Rcout << tmp_global_stats << std::endl;
    // Rcout << steps.at(i) << std::endl;
    
    // bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    // if(is_converged){
    //   break;
    // }
    
  }
  Rcout << "Successfully found the starting value" << std::endl;
  Rcout << "Start MLE" << std::endl; 
  return(coefs);
}

// [[Rcpp::export]]     
std::vector<arma::vec> t_mle_estimation_p(std::vector<arma::mat> networks,
                                                 std::vector<std::string> terms, 
                                                 int n_actors,
                                                 int n_proposals,
                                                 int n_proposals_burn_in,
                                                 int seed,
                                                 int number_networks, 
                                                 std::vector<std::vector<arma::mat>> &data_lists, 
                                                 std::vector<int> &type_list,
                                                 arma::vec &beg_coef,
                                                 bool mh,
                                                 List cluster,
                                                 List data_lists_par,
                                                 List networks_par,
                                                 int max_it = 10, 
                                                 double tol = 0.001, 
                                                 bool start_with_empty_net = true) {
  
  // Set up objects
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_pos;
  std::unordered_map< int, std::set<int>> pos_edges;
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_neg = all_edges_pos;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    all_edges_pos.push_back(mat_to_map(networks.at(i),1, n_actors));
    all_edges_neg.push_back(mat_to_map(networks.at(i),-1, n_actors));
  }
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean(1,terms.size()), tmp_var(terms.size(),terms.size());
  arma::vec update_coef,tmp_coef; 
  std::vector<arma::mat> all_stats;
  arma::mat trying,tmp_stats(number_networks,terms.size()), sum_stats(number_networks,terms.size());
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::mat global_stats_per_network(data_lists.size(),type_list.size());
  arma::vec sum_global_stats(type_list.size()), tmp(type_list.size());
  // Rcout << all_edges_pos.size() << std::endl;
  // Rcout << all_edges_neg.size() << std::endl;
  
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    tmp = count_global_statistic(all_edges_pos.at(i),
                                 all_edges_neg.at(i),
                                 n_actors, 
                                 data_lists.at(i), 
                                 type_list,
                                 functions) ;
    // Rcout << tmp << std::endl;
    // Rcout << at_zero << std::endl;
    // Rcout << tmp + at_zero << std::endl;
    
    tmp += at_zero;
    global_stats_per_network.row(i) = tmp.as_row();
    sum_global_stats += tmp;
  }
  coefs.push_back(beg_coef);
  Function parallel_sample("parallel_sample");
  List res;
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl;
    tmp_coef = coefs.at(i);
    // for each network simulate the statistics and estimate their mean and covariance 
    res = parallel_sample(tmp_coef,terms, n_actors,
                          networks_par,data_lists_par,type_list,mh,
                          n_proposals,n_proposals_burn_in,
                          seed +i,number_networks,cluster,start_with_empty_net);
    // Rcout << n.at(0) << std::endl;
    NumericVector tmp_mean_1 = res.at(0);
     
    NumericMatrix tmp_var_1 = res.at(1);
    NumericMatrix sum_global_stats_1 = res.at(2);
    tmp_mean = arma::mat(tmp_mean_1.begin(),  1, tmp_mean_1.size(),false);
    tmp_var = arma::mat(tmp_var_1.begin(), tmp_var_1.nrow(), tmp_var_1.ncol(), false);
    sum_stats = arma::mat(sum_global_stats_1.begin(), sum_global_stats_1.nrow(), sum_global_stats_1.ncol(), false);
    update_coef = tmp_coef + arma::solve(tmp_var,sum_global_stats-tmp_mean.t());
    coefs.push_back(update_coef);
    
    // Set all values back to zero 
    tmp_var.zeros();
    tmp_mean.zeros();
    sum_stats.zeros();
    
    
    // Rcout << tmp_mean << std::endl;
    // Rcout << tmp_global_stats << std::endl;
    // Rcout << steps.at(i) << std::endl;
    
    bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    if(is_converged){
      break;
    }
    
  }
  return(coefs);
}


// [[Rcpp::export]]     
std::vector<arma::vec> t_mle_estimation(std::vector<arma::mat> networks,
                                        std::vector<std::string> terms, 
                                        int n_actors,
                                        int n_proposals,
                                        int n_proposals_burn_in,
                                        int seed,
                                        int number_networks, 
                                        std::vector<std::vector<arma::mat>> &data_lists, 
                                        std::vector<int> &type_list,
                                        arma::vec &beg_coef,
                                        bool mh,
                                        int max_it = 10, 
                                        double tol = 0.001, 
                                        bool start_with_empty_net = true) {
  
  // Set up objects
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_pos;
  std::unordered_map< int, std::set<int>> pos_edges;
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_neg = all_edges_pos;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    all_edges_pos.push_back(mat_to_map(networks.at(i),1, n_actors));
    all_edges_neg.push_back(mat_to_map(networks.at(i),-1, n_actors));
  }
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean(1,terms.size()), tmp_var(terms.size(),terms.size());
  arma::vec update_coef,tmp_coef; 
  std::vector<arma::mat> all_stats;
  arma::mat trying,tmp_stats(number_networks,terms.size()), sum_stats(number_networks,terms.size());
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::mat global_stats_per_network(data_lists.size(),type_list.size());
  arma::vec sum_global_stats(type_list.size()), tmp(type_list.size());
  // Rcout << all_edges_pos.size() << std::endl;
  // Rcout << all_edges_neg.size() << std::endl;
  
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    tmp = count_global_statistic(all_edges_pos.at(i),
                                 all_edges_neg.at(i),
                                 n_actors, 
                                 data_lists.at(i), 
                                 type_list,
                                 functions) ;
    // Rcout << tmp << std::endl;
    // Rcout << at_zero << std::endl;
    // Rcout << tmp + at_zero << std::endl;
    
    tmp += at_zero;
    global_stats_per_network.row(i) = tmp.as_row();
    sum_global_stats += tmp;
  }
  coefs.push_back(beg_coef);
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl;
    tmp_coef = coefs.at(i);
    // for each network simulate the statistics and estimate their mean and covariance 
    for(unsigned int j = 0; j< data_lists.size();j +=1){
      if(start_with_empty_net){
        // a) Sample
        tmp_stats = simulate_networks_stats(tmp_coef,terms, n_actors,
                                            data_lists.at(j), type_list,mh,
                                            n_proposals,n_proposals_burn_in,
                                            seed +i,number_networks);
        all_stats.push_back(tmp_stats);
        sum_stats += tmp_stats;
        // b) Estimate mean 
        tmp_mean += arma::mean(tmp_stats,0);
        // c) Estimate variance 
        tmp_var += arma::cov(tmp_stats);
      } else {
        tmp_stats=simulate_networks_intern(tmp_coef,
                                           all_edges_pos.at(j),
                                           all_edges_neg.at(j),
                                           global_stats_per_network.row(j).t(),
                                           mh,
                                           functions,
                                           n_actors,
                                           data_lists.at(j),
                                           type_list,
                                           n_proposals,
                                           n_proposals_burn_in,
                                           seed +i,
                                           number_networks);
        all_stats.push_back(tmp_stats);
        sum_stats += tmp_stats;
        // b) Estimate mean 
        tmp_mean += arma::mean(tmp_stats,0);
        // c) Estimate variance 
        tmp_var += arma::cov(tmp_stats);
      }
    }
    tmp_global_stats = sum_global_stats;
    update_coef = tmp_coef + arma::solve(tmp_var,tmp_global_stats-tmp_mean.t());
    coefs.push_back(update_coef);
    
    // Set all values back to zero 
    tmp_var.zeros();
    tmp_mean.zeros();
    sum_stats.zeros();
    
    
    // Rcout << tmp_mean << std::endl;
    // Rcout << tmp_global_stats << std::endl;
    // Rcout << steps.at(i) << std::endl;
    
    bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    if(is_converged){
      break;
    }
    
  }
  return(coefs);
}


// [[Rcpp::export]]     
List t_est_var_p(std::vector<arma::mat> networks,
               std::vector<std::string> terms, 
               int n_actors,
               int n_proposals,
               int n_proposals_burn_in,
               int seed,
               int number_networks, 
               std::vector<std::vector<arma::mat>> &data_lists, 
               std::vector<int> &type_list,
               arma::vec &coef,
               bool mh,
               List cluster,
               List data_lists_par,
               List networks_par,
               bool start_with_empty_net = true) {
  
  // Set up objects
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_pos;
  std::unordered_map< int, std::set<int>> pos_edges;
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_neg = all_edges_pos;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    all_edges_pos.push_back(mat_to_map(networks.at(i),1, n_actors));
    all_edges_neg.push_back(mat_to_map(networks.at(i),-1, n_actors));
  } 
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean(1,terms.size()), tmp_var(terms.size(),terms.size());
  arma::vec update_coef,tmp_coef; 
  // std::vector<arma::mat> all_stats;
  arma::mat trying,tmp_stats(number_networks,terms.size()), sum_stats(number_networks,terms.size()), sum_stats_raw(number_networks,terms.size());
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::mat global_stats_per_network(data_lists.size(),type_list.size());
  arma::vec sum_global_stats(type_list.size()), tmp(type_list.size());
  

  
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    tmp = count_global_statistic(all_edges_pos.at(i),
                                 all_edges_neg.at(i),
                                 n_actors, 
                                 data_lists.at(i), 
                                 type_list,
                                 functions) ;
    // Rcout << tmp << std::endl;
    // Rcout << at_zero << std::endl;
    // Rcout << tmp + at_zero << std::endl;
    
    tmp += at_zero;
    global_stats_per_network.row(i) = tmp.as_row();
    sum_global_stats += tmp;
  } 
  
  Function parallel_sample("parallel_sample");
  List res;
  res =  parallel_sample(coef,terms, n_actors,
                         networks_par,data_lists_par,type_list,mh,
                         n_proposals,n_proposals_burn_in,
                         seed,number_networks,cluster,start_with_empty_net);
  // Rcout << n.at(0) << std::endl;
  NumericVector tmp_mean_1 = res.at(0);
   
  NumericMatrix tmp_var_1 = res.at(1);
  NumericMatrix sum_global_stats_1 = res.at(2);
  tmp_mean = arma::mat(tmp_mean_1.begin(),  1, tmp_mean_1.size(),false);
  tmp_var = arma::mat(tmp_var_1.begin(), tmp_var_1.nrow(), tmp_var_1.ncol(), false);
  sum_stats = arma::mat(sum_global_stats_1.begin(), sum_global_stats_1.nrow(), sum_global_stats_1.ncol(), false);
  
  
  // Get the std error of the sufficient statistics 
  arma::mat std_err = sqrt(diagvec(arma::cov(sum_stats)));
  // Rcout << std_err << std::endl;
  
  // Get the means of the sufficient statistics 
  arma::mat t_vals = (tmp_mean - sum_global_stats.t())/(std_err.t());
  sum_stats_raw = sum_stats;
  // Normalize the sampled statistics 
  sum_stats.each_row() -= sum_global_stats.t();
  return(List::create(Named("var") = arma::inv(arma::cov(sum_stats)),
                      // _["var_alt"] = arma::inv(arma::cov(sum_stats)),
                      // _["var_alt_2"] = arma::inv(arma::cov(tmp_stats)),
                      _["stats_mean"] = sum_stats, 
                      _["stats_raw"] = sum_stats_raw, 
                      _["obs"] = sum_global_stats,
                      _["t_vals"] = t_vals));
  
} 

// [[Rcpp::export]]     
List t_est_var(std::vector<arma::mat> networks,
                                        std::vector<std::string> terms, 
                                        int n_actors,
                                        int n_proposals,
                                        int n_proposals_burn_in,
                                        int seed,
                                        int number_networks, 
                                        std::vector<std::vector<arma::mat>> &data_lists, 
                                        std::vector<int> &type_list,
                                        arma::vec &coef,
                                        bool mh,
                                        bool start_with_empty_net = true) {
  
  // Set up objects
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_pos;
  std::unordered_map< int, std::set<int>> pos_edges;
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_neg = all_edges_pos;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    all_edges_pos.push_back(mat_to_map(networks.at(i),1, n_actors));
    all_edges_neg.push_back(mat_to_map(networks.at(i),-1, n_actors));
  }
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean(1,terms.size()), tmp_var(terms.size(),terms.size());
  arma::vec update_coef,tmp_coef; 
  // std::vector<arma::mat> all_stats;
  arma::mat trying,tmp_stats(number_networks,terms.size()), sum_stats(number_networks,terms.size()),sum_stats_raw(number_networks,terms.size());
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::mat global_stats_per_network(data_lists.size(),type_list.size());
  arma::vec sum_global_stats(type_list.size()), tmp(type_list.size());

  for(unsigned int i = 0; i< data_lists.size();i +=1){
    tmp = count_global_statistic(all_edges_pos.at(i),
                                 all_edges_neg.at(i),
                                 n_actors, 
                                 data_lists.at(i), 
                                 type_list,
                                 functions) ;
    // Rcout << tmp << std::endl;
    // Rcout << at_zero << std::endl;
    // Rcout << tmp + at_zero << std::endl;
    
    tmp += at_zero;
    global_stats_per_network.row(i) = tmp.as_row();
    sum_global_stats += tmp;
  }
  
  // for each network simulate the statistics and estimate their mean and covariance 
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    if(start_with_empty_net){
      // a) Sample
      tmp_stats = simulate_networks_stats(coef,terms, n_actors,
                                          data_lists.at(i), type_list,mh,
                                          n_proposals,n_proposals_burn_in,
                                          seed +i,number_networks);
      sum_stats += tmp_stats;
      // b) Estimate mean 
      tmp_mean += arma::mean(tmp_stats,0);
      // c) Estimate variance 
      tmp_var += arma::cov(tmp_stats);
    } else {
      tmp_stats=simulate_networks_intern(coef,
                                         all_edges_pos.at(i),
                                         all_edges_neg.at(i),
                                         global_stats_per_network.row(i).t(),
                                         mh,
                                         functions,
                                         n_actors,
                                         data_lists.at(i),
                                         type_list,
                                         n_proposals,
                                         n_proposals_burn_in,
                                         seed +i,
                                         number_networks);
      sum_stats += tmp_stats;
      // b) Estimate mean 
      tmp_mean += arma::mean(tmp_stats,0);
      // c) Estimate variance 
      tmp_var += arma::cov(tmp_stats);
    }
  }
  // Get the std error of the sufficient statistics 
  arma::mat std_err = sqrt(diagvec(arma::cov(sum_stats)));
  // Rcout << std_err << std::endl;
  
  // Get the means of the sufficient statistics 
  arma::mat t_vals = (tmp_mean - sum_global_stats.t())/(std_err.t());
  
  // Rcout << tmp_var << std::endl;
  List res; 
  sum_stats_raw = sum_stats;
  // Normalize the sampled statistics 
  sum_stats.each_row() -= sum_global_stats.t();
  return(List::create(Named("var") = arma::inv(arma::cov(sum_stats)),
                      // _["var_alt"] = arma::inv(arma::cov(sum_stats)),
                      // _["var_alt_2"] = arma::inv(arma::cov(tmp_stats)),
                      _["stats_mean"] = sum_stats, 
                      _["stats_raw"] = sum_stats_raw, 
                      _["obs"] = sum_global_stats,
                      _["t_vals"] = t_vals));
  
}

// [[Rcpp::export]]     
std::vector<arma::vec> t_rm_mle_estimation(std::vector<arma::mat> networks,
                                        std::vector<std::string> terms, 
                                        int n_actors,
                                        int n_proposals,
                                        int n_proposals_burn_in,
                                        int seed,
                                        int number_networks,
                                        std::vector<std::vector<arma::mat>> &data_lists, 
                                        std::vector<int> &type_list,
                                        arma::vec &beg_coef,
                                        bool mh,
                                        arma::mat D,
                                        double c, 
                                        int max_it = 10, 
                                        double tol = 0.001, 
                                        bool start_with_empty_net = true) {
  
  // Set up objects
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_pos;
  std::unordered_map< int, std::set<int>> pos_edges;
  std::vector<std::unordered_map< int, std::set<int>>> all_edges_neg = all_edges_pos;
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    all_edges_pos.push_back(mat_to_map(networks.at(i),1, n_actors));
    all_edges_neg.push_back(mat_to_map(networks.at(i),-1, n_actors));
  } 
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean(1,terms.size()), tmp_var(terms.size(),terms.size());
  arma::vec update_coef,tmp_coef; 
  std::vector<arma::mat> all_stats;
  arma::mat trying,tmp_stats(number_networks,terms.size()), sum_stats(number_networks,terms.size());
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  arma::vec tmp_global_stats;
  arma::mat global_stats_per_network(data_lists.size(),type_list.size());
  arma::vec sum_global_stats(type_list.size()), tmp(type_list.size());
  
  for(unsigned int i = 0; i< data_lists.size();i +=1){
    tmp = count_global_statistic(all_edges_pos.at(i),
                                 all_edges_neg.at(i),
                                 n_actors, 
                                 data_lists.at(i), 
                                 type_list,
                                 functions) ;
    
    tmp += at_zero;
    global_stats_per_network.row(i) = tmp.as_row();
    sum_global_stats += tmp;
  } 
  coefs.push_back(beg_coef);
  
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl;
    tmp_coef = coefs.at(i);
    // for each network simulate the statistics and estimate their mean and covariance 
    for(unsigned int j = 0; j< data_lists.size();j ++){
      if(start_with_empty_net){
        // a) Sample
        tmp_stats = simulate_networks_stats(tmp_coef,terms, n_actors,
                                            data_lists.at(j), type_list,mh,
                                            n_proposals,n_proposals_burn_in,
                                            seed +i,number_networks);
        all_stats.push_back(tmp_stats);
        sum_stats += tmp_stats;
        // b) Estimate mean 
        tmp_mean += arma::mean(tmp_stats,0);
        // c) Estimate variance 
        tmp_var += arma::cov(tmp_stats);
      } else { 
        tmp_stats=simulate_networks_intern(tmp_coef,
                                           all_edges_pos.at(j),
                                           all_edges_neg.at(j),
                                           global_stats_per_network.row(j).t(),
                                           mh,
                                           functions,
                                           n_actors,
                                           data_lists.at(j),
                                           type_list,
                                           n_proposals,
                                           n_proposals_burn_in,
                                           seed +i,
                                           number_networks);
        all_stats.push_back(tmp_stats);
        sum_stats += tmp_stats;
        // b) Estimate mean 
        tmp_mean += arma::mean(tmp_stats,0);
        // c) Estimate variance 
        tmp_var += arma::cov(tmp_stats);
      } 
    }
    update_coef = tmp_coef +  std::pow((i+1),-c)*D*(sum_global_stats-tmp_mean.t());
    coefs.push_back(update_coef);
    
    // Set all values back to zero 
    tmp_var.zeros();
    tmp_mean.zeros();
    sum_stats.zeros();
    
    bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    if(is_converged){
      break;
    } 
    
  } 
  return(coefs);
} 


// [[Rcpp::export]]
List est_var(arma::mat network ,
             std::vector<std::string> terms, 
             int n_actors,
             int n_proposals,
             int n_proposals_burn_in,
             int seed,
             int number_networks, 
             std::vector<arma::mat> &data_list, 
             std::vector<int> &type_list,
             arma::vec coef,
             bool mh,
             bool start_with_empty_net = true) {
  
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  arma::mat tmp_stats;

  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg,
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  Rcout << "Start" << std::endl;
  
  if(start_with_empty_net){
    tmp_stats = simulate_networks_stats(coef,terms, n_actors,
                                        data_list, type_list,mh,
                                        n_proposals,n_proposals_burn_in, 
                                        seed,number_networks);
  } else {
    tmp_stats = simulate_networks_intern(coef,edges_pos, edges_neg,global_stats,mh, functions, n_actors,
                                         data_list, type_list,
                                         n_proposals,n_proposals_burn_in,
                                         seed,number_networks);
  }
  Rcout << "Done" << std::endl;
  // Rcout << global_stats << std::endl;
  
  arma::mat tmp_var = arma::cov(tmp_stats);
  // Get the std error of the sufficient statistics 
  arma::mat std_err = sqrt(diagvec(tmp_var));
  // Rcout << std_err << std::endl;
  
  // Get the means of the sufficient statistics 
  arma::mat mean = arma::mean(tmp_stats,0);
  // Rcout << mean << std::endl;
  arma::mat t_vals = (mean - global_stats.t())/(std_err.t());
  
  // Rcout << tmp_var << std::endl;
  List res; 
  // Normalize the sampled statistics 
  tmp_stats.each_row() -= global_stats.t();
  return(List::create(Named("var") = arma::inv(tmp_var),
                      _["stats"] = tmp_stats, 
                      _["t_vals"] = t_vals));
}


// Here the mle estimation is done  according to Snijders    
// [[Rcpp::export]]     
std::vector<arma::vec> rm_mle_estimation(arma::mat network ,
                                         std::vector<std::string> terms, 
                                         int n_actors,
                                         int n_proposals,
                                         int n_proposals_burn_in,
                                         int seed,
                                         int number_networks, 
                                         std::vector<arma::mat> &data_list, 
                                         std::vector<int> &type_list,
                                         arma::vec &beg_coef,
                                         bool mh,
                                         arma::mat D,
                                         double c, 
                                         int max_it = 10, 
                                         double tol = 0.001, 
                                         bool start_with_empty_net = true) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  std::vector<arma::vec> coefs;
  arma::mat tmp_mean,tmp_stats, tmp_var;
  arma::vec update_coef,tmp_coef; 
  
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);  
  List tmp_samples;
  
  arma::vec at_zero;
  at_zero = eval_at_empty_network(terms, n_actors);
  
  arma::vec global_stats(count_global_statistic(edges_pos,
                                                edges_neg,
                                                n_actors, 
                                                data_list, 
                                                type_list,
                                                functions));
  global_stats = global_stats + at_zero;
  // Set \theta_0
  coefs.push_back(beg_coef);
  Rcout << "Phase 1 finished" << std::endl;
  for(int i = 0; i <(max_it);i ++) {
    Rcout << "Iteration:" + std::to_string(i+1) << std::endl; 
    tmp_coef = coefs.at(i);
    if(start_with_empty_net){
      tmp_stats = simulate_networks_stats(tmp_coef,terms, n_actors,
                                          data_list, type_list,mh,
                                          n_proposals,n_proposals_burn_in, 
                                          seed +i,number_networks);
    } else { 
      tmp_stats = simulate_networks_intern(tmp_coef,
                                           edges_pos, 
                                           edges_neg,
                                           global_stats,
                                           mh, 
                                           functions, 
                                           n_actors,
                                           data_list, 
                                           type_list,
                                           n_proposals,
                                           n_proposals_burn_in,
                                           seed +i,
                                           number_networks);
    }
    
    // Rcout << tmp_stats << std::endl;
    // Rcout << tmp_coef << std::endl;
    // tmp_stats =as<arma::mat> (tmp_samples["stats"]);
    // Rcout << tmp_stats << std::endl;
    // Rcout << tmp_var << std::endl;
    
    // Rcout << std::pow((i+1),-c) << std::endl;
    tmp_mean = arma::mean(tmp_stats,0);
    
    update_coef = tmp_coef +  std::pow((i+1),-c)*D*(global_stats-tmp_mean.t());
    
    // update_coef = tmp_coef + arma::inv(tmp_var)*(global_stats-tmp_mean.t());
    coefs.push_back(update_coef);
    // Rcout << tmp_mean << std::endl;
    // Rcout << update_coef << std::endl;
    
    bool is_converged = std::sqrt(arma::accu(arma::square(update_coef - tmp_coef))) < tol;
    if(is_converged){
      break;
    }
    // coefs.push_back(arma::var(tmp_stats));
  }
  return(coefs);
}





// [[Rcpp::export]]
List preprocess_pseudo_lh_new(arma::mat network, 
                              int n_actors,
                              std::vector<std::string> terms,
                              std::vector<arma::mat> &data_list, 
                              std::vector<int> &type_list) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  arma::mat change_stat;
  arma::vec tmp_row;
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  
  // Set up objects
  arma::mat tmp_mat(n_actors, n_actors);
  arma::vec tmp_vec(n_actors);
  tmp_vec.fill(R_NaN);
  //For now we only have two entries in the lists 
  std::vector<arma::mat> pos_change;
  std::vector<arma::mat> neg_change;
  std::vector<arma::mat> zero_change;
  arma::mat  mat_tmp(n_actors, n_actors);
  mat_tmp.diag() = tmp_vec;
  mat_tmp(arma::trimatl_ind(size(mat_tmp))).fill(R_NaN);
  // Rcout << mat_tmp << std::endl;
  
  for(unsigned int i = 0; i <(terms.size()); i = i + 1 ) {
    pos_change.push_back(mat_tmp);
    neg_change.push_back(mat_tmp);
    zero_change.push_back(mat_tmp);
  }
  
  arma::umat pos_response(n_actors, n_actors);
  arma::umat neg_response(n_actors, n_actors);
  arma::umat zero_response(n_actors, n_actors);
  arma::vec global_statistics(terms.size());
  
  int present_val;
  
  // Fill the response matrices
  arma::mat ind_network(n_actors,n_actors);
  ind_network.fill(1);
  pos_response = network == ind_network;
  ind_network.fill(-1);
  neg_response = network == ind_network;
  ind_network.fill(0);
  zero_response = network == ind_network;
  // Generate change statistic function from the terms 
  std::vector<arma::vec(*)( const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &present_val,
                        int &proposed_change, 
                        arma::mat &data, 
                        int &type)> functions;
  functions = change_statistics_generate(terms);
  
  // For each entry in network we need to calculate the change stats
  // (in all directions)
  for(int i = 1; i <=(n_actors); i = i + 1 ) {
    for(int j = i +1; j <=(n_actors); j = j + 1 ) {
      if(i == j){
        continue;
      }
      // Get the present value of i,j
      if(edges_pos.at(i).count(j)){
        present_val = 1;
      } else if(edges_neg.at(i).count(j)){
        present_val = -1;
      } else {
        present_val = 0;
      }
      // Calculate the change statstics 
      
      
      change_stat =   calculate_change_stats(i,
                                             j,
                                             n_actors,
                                             edges_pos,
                                             edges_neg,  
                                             data_list, 
                                             type_list,
                                             functions);
      // Rcout << change_stat << std::endl;
      
      // We now need to calculate delta(y_ij^+), delta(y_ij^-), delta(y_ij^0)
      if(present_val == 0){
        for(unsigned int n = 0; n <(functions.size()); n = n + 1 ) {
          pos_change.at(n).at(i-1,j-1) = change_stat(n,0);
          neg_change.at(n).at(i-1,j-1) = change_stat(n,1);
          zero_change.at(n).at(i-1,j-1) = 0;
        }
      }
      // Rcout << "1" << std::endl;
      
      if(present_val == 1){
        for(unsigned int n = 0; n <(functions.size()); n = n + 1 ) {
          pos_change.at(n).at(i-1,j-1) = 0;
          neg_change.at(n).at(i-1,j-1) = -change_stat(n,0);
          zero_change.at(n).at(i-1,j-1) = change_stat(n,2);
        }
      }
      // Rcout << "1" << std::endl;
      if(present_val == -1){
        for(unsigned int n = 0; n <(functions.size()); n = n + 1 ) {
          pos_change.at(n).at(i-1,j-1) = -change_stat(n,2);
          neg_change.at(n).at(i-1,j-1) = 0;
          zero_change.at(n).at(i-1,j-1) = -change_stat(n,1);
        }
      }
      // Rcout << "1" << std::endl;
    } 
  }
  return(List::create(Named("pos_changes") = pos_change,
                      _["neg_changes"] = neg_change,
                      _["zero_changes"] = zero_change,
                      _["pos_response"] = pos_response,
                      _["neg_response"] =neg_response,
                      _["zero_response"] =zero_response));
}

//' Count the number of dyadwise shared partners
//' 
//' This function returns for each dyad the number of shared positive partners.
//'
// @param network A matrix of the adjacency matrix of a signed network. 
// @param n_actors A numeric value of the number of actors in the network. 
// @export
// [[Rcpp::export]]
arma::vec count_dyadwise_shared_partner_pos(arma::mat network, 
                              int n_actors) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  arma::mat change_stat;
  arma::vec res(n_actors*(n_actors-1)/2);
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  std::set<int> intersection;
  int n = 0;
  for(int i = 1; i <=(n_actors); i = i + 1 ) {
    for(int j = i +1; j <=(n_actors); j = j + 1 ) {
      if(i == j){
        continue;
        Rcout << "1" << std::endl;
      }
      // Step 1: Get all connections of tmp_random_i
      std::set<int> connections_of_i = edges_pos.at(i);
      // Step 2: Get all connections of tmp_random_j
      std::set<int> connections_of_j = edges_pos.at(j);
      intersection.clear();    
      std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
                            std::begin(connections_of_j), std::end(connections_of_j),
                            std::inserter(intersection, std::begin(intersection)));
      res.at(n) = intersection.size();
      n +=1;
    } 
  }
  return(res);
}
//' Count the number of positive edgewise shared partners
//' 
//' This function returns for each dyad the number of shared positive partners.
//'
//' @param network A matrix of the adjacency matrix of a signed network. 
//' @param n_actors A numeric value of the number of actors in the network. 
//' @param n_edges A numeric value of the number of edges in the network. 
// @export
// [[Rcpp::export]]
arma::vec count_edgewise_shared_partner_pos(arma::mat network, 
                                   int n_actors, 
                                   int n_edges) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  // std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  arma::mat change_stat;
  arma::vec res(n_edges);
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  // edges_neg = mat_to_map(network,-1, n_actors);
  std::set<int> intersection;
  int n = 0;
  for(int i = 1; i <=(n_actors); i = i + 1 ) {
    for(int j = i +1; j <=(n_actors); j = j + 1 ) {
      if(edges_pos.at(i).count(j)){
        // Step 1: Get all connections of tmp_random_i
        std::set<int> connections_of_i = edges_pos.at(i);
        // Step 2: Get all connections of tmp_random_j
        std::set<int> connections_of_j = edges_pos.at(j);
        intersection.clear();    
        std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
                              std::begin(connections_of_j), std::end(connections_of_j),
                              std::inserter(intersection, std::begin(intersection)));
        res.at(n) = intersection.size();
        n +=1;
      }
    } 
  }
  return(res);
}
//' Count the number of edgewise shared enemies
//' 
//' This function returns for each edge the number of shared enemies.
//'
//' @param network A matrix of the adjacency matrix of a signed network. 
//' @param n_actors A numeric value of the number of actors in the network. 
//' @param n_edges A numeric value of the number of edges in the network. 
// @export
// [[Rcpp::export]]
arma::vec count_edgewise_shared_enemies(arma::mat network, 
                                            int n_actors, 
                                            int n_edges) {
  // Set up objects
  std::unordered_map< int, std::set<int>> edges_pos;
  std::unordered_map< int, std::set<int>> edges_neg = edges_pos;
  arma::mat change_stat;
  arma::vec res(n_edges);
  // Convert the matrix to two unordered_map objects 
  edges_pos = mat_to_map(network,1, n_actors);
  edges_neg = mat_to_map(network,-1, n_actors);
  std::set<int> intersection;
  int n = 0;
  for(int i = 1; i <=(n_actors); i = i + 1 ) {
    for(int j = i +1; j <=(n_actors); j = j + 1 ) {
      if(edges_pos.at(i).count(j)){
        // Step 1: Get all connections of tmp_random_i
        std::set<int> connections_of_i = edges_neg.at(i);
        // Step 2: Get all connections of tmp_random_j
        std::set<int> connections_of_j = edges_neg.at(j);
        intersection.clear();    
        std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
                              std::begin(connections_of_j), std::end(connections_of_j),
                              std::inserter(intersection, std::begin(intersection)));
        res.at(n) = intersection.size();
        n +=1;
      }
    } 
  } 
  return(res);
} 
// 
// // [[Rcpp::export]]
// arma::vec trying_out_2(  std::vector<std::string> terms,arma::mat network, 
//                          std::vector<List> &data_list, 
//                          std::vector<int> &type_list,
//                          int n_actors = 10) {
//   std::vector<arma::vec(*)( std::unordered_map< int, std::set<int>> &edges_pos,
//                         std::unordered_map< int, std::set<int>> &edges_neg,
//                         int &n_actors,
//                         int &present_val,
//                         int &proposed_change, 
//                         List &data, 
//                         int &type)> functions;
//   functions = change_statistics_generate(terms);
//   
//   arma::mat res(n_actors,n_actors);
//   res.fill(0);
//   std::unordered_map< int, std::set<int>> res_plus;
//   res_plus = mat_to_map(network,1,n_actors);
//   std::unordered_map< int, std::set<int>> res_minus;
//   res_minus = mat_to_map(network,-1,n_actors);
//   res = count_global_statistic(res_plus,
//                                res_minus,
//                                n_actors, 
//                                data_list, 
//                                type_list,
//                                functions) ;
//   
//   return(res);
// }
// 
