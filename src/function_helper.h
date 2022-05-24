// Defines a header file containing function signatures for functions in src/

// Protect signatures using an inclusion guard.
#ifndef function_helper_H
#define function_helper_H

std::unordered_map< int, std::set<int>> mat_to_map(arma::mat mat, int type = 1, int n_actors = 10);
arma::mat map_to_mat(std::unordered_map< int, std::set<int>> &adj_plus,
                     std::unordered_map< int, std::set<int>> &adj_minus,
                     int n_actors = 10);
arma::mat trying_out(arma::mat network,
                     int n_actors = 10);
void set_seed(double seed); 
arma::vec simulate_numbers(int min = 0, int max = 10, int number = 100) ;
#endif
