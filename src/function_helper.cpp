// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <random>
#include <set>
#include <unordered_map>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// This is a function to transform one matrix to one maps (either for positive or negative ties)
// [[Rcpp::export]]
std::unordered_map< int, std::set<int>> mat_to_map(arma::mat mat, int type = 1, int n_actors = 10) {
  std::unordered_map< int, std::set<int>> adj_e;
  for (int i = 1; i <= n_actors; i++){
    adj_e[i] = std::set<int>();
  }
  arma::vec tmp_row;
  for (int i = 1; i <= n_actors; i++){
    tmp_row = mat.col(i-1);
    arma::uvec ids = find(tmp_row == type) + 1; // Find indices
    adj_e[i] = std::set<int>(ids.begin(),ids.end());
  }
  return(adj_e);
}

// This is a function to transform two maps (one for positive and one for negative ties) to one matrix
arma::mat map_to_mat(std::unordered_map< int, std::set<int>> &adj_plus,
                     std::unordered_map< int, std::set<int>> &adj_minus,
                     int n_actors = 10) {
  arma::mat res(n_actors,n_actors);
  res.fill(0);
  arma::mat tmp_row(n_actors,1);
  NumericVector tmp_vec_plus(n_actors);
  NumericVector tmp_vec_minus(n_actors);
  for (int i = 1; i < n_actors; i++){
    // Set the vector to 0
    tmp_row.fill(0); 
    tmp_vec_plus = adj_plus[i];
    tmp_vec_minus = adj_minus[i];
    for (int j = 0; j < tmp_vec_plus.length() ; j++){
      // tmp_row[tmp_vec_plus(j)-1] = 1;
      tmp_row.at(tmp_vec_plus.at(j)-1) = 1;
    }
    for (int j = 0; j < tmp_vec_minus.length(); j++){
      // tmp_row[tmp_vec_minus(j)-1] = -1;
      tmp_row.at(tmp_vec_minus.at(j)-1) = - 1;
    }
    res.col(i-1) = tmp_row.as_col();
    // if(i == 2){
    //   Rcout << "Mean " <<  tmp_row.head(10) << std::endl;
    //   // Rcout << "Mean " <<  tmp_row.head(20) << std::endl;
    //   Rcout << "tmp_vec_plus " <<  tmp_vec_plus << std::endl;
    //   Rcout << "tmp_vec_plus " <<  tmp_vec_minus << std::endl;
    //   Rcout << "Max " <<  tmp_row.n_elem << std::endl;
    // }
  }
  return(res);
}
// [[Rcpp::export]]
arma::mat trying_out(arma::mat network,
                     int n_actors = 10) {
  arma::mat res(n_actors,n_actors);
  res.fill(0);
  std::unordered_map< int, std::set<int>> res_plus;
  res_plus = mat_to_map(network,1,n_actors);
  std::unordered_map< int, std::set<int>> res_minus;
  res_minus = mat_to_map(network,-1,n_actors);
  res = map_to_mat(res_plus, res_minus,n_actors);
  return(res);
}

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]
arma::vec simulate_numbers(int min = 0, int max = 10, int number = 100) {
  int range = max - min +1;
  arma::vec res(number);
  res.randu();
  res = arma::trunc(res*range) + min;
  return(res);
}
