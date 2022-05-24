// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <random>
#include <set>
#include <unordered_map>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp; 

// The dimension of change_stat equals the number of terms in the model
// Finished Terms: 
// 0: Number of positive edges 
// 1: Number of negative edges 
// Terms in Planning:
// 2: Positive Balanced triangles: i-k-j-i where all connections are positive 
// 3: Negative Balanced triangles: i-k-j-i where i-k and k-j negative and i-j positive 
// 4: k-Positive Balanced triangles: k common partners between i and j where all connections are positive 
// 5: k-Negative Balanced triangles: k common partners between i and j where only i-j is positive
// 6: GW 5 and 6 
// 7. 4-7 for negative and positive two-stars 
// 8. 4-7 for negative and positive degree 
// 9. Mean correlation of positive and negative degrees
// 10. Number of isolates (can be in general or specific to negative or positive networks)
// 11. Covariate effects same as in the standard ergm but same conditions as for 11. apply (Only dyadic)
// The dimension needs to be set according to the dimension of the change stats  

//  Calculate and return  
//  delta(present_val)^(0->1) 
//  delta(present_val)^(0->-1) 
//  delta(present_val)^(1->-1) 
// Change statistics of edges_pos () 
// 0: Number of positive edges 
arma::vec stat_edges_pos( const std::unordered_map< int, std::set<int>> &edges_pos,
                          const std::unordered_map< int, std::set<int>> &edges_neg,
                          int &n_actors,
                          int &tmp_random_i,
                          int &tmp_random_j, 
                          arma::mat &data, 
                          int &type){
  arma::vec res(3);
  res = {1,0,-1};
  return(res);
}

arma::vec stat_edges( const std::unordered_map< int, std::set<int>> &edges_pos,
                      const std::unordered_map< int, std::set<int>> &edges_neg,
                      int &n_actors,
                      int &tmp_random_i,
                      int &tmp_random_j, 
                      arma::mat &data, 
                      int &type){
  arma::vec res(3);
  res = {1,1,0};
  return(res);
}
// Number of negative edges 
arma::vec stat_edges_neg(  const std::unordered_map< int, std::set<int>> &edges_pos,
                           const std::unordered_map< int, std::set<int>> &edges_neg,
                           int &n_actors,
                           int &tmp_random_i,
                           int &tmp_random_j, 
                           arma::mat &data, 
                           int &type){
  arma::vec res(3);
  res = {0,1,1};
  return(res);
}
arma::vec stat_isolates( const std::unordered_map< int, std::set<int>> &edges_pos,
                          const std::unordered_map< int, std::set<int>> &edges_neg,
                          int &n_actors,
                          int &tmp_random_i,
                          int &tmp_random_j, 
                          arma::mat &data, 
                          int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size() + edges_neg.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size() + edges_neg.at(tmp_random_j).size();
  // If the edge is already there, we need to substract one of the degrees
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
 
  res.at(0) = -(degree_i == 0) -(degree_j == 0);
  res.at(1) = res.at(0) ;
  res.at(2) = 0;
  return(res);
}

arma::vec stat_isolates_pos( const std::unordered_map< int, std::set<int>> &edges_pos,
                         const std::unordered_map< int, std::set<int>> &edges_neg,
                         int &n_actors,
                         int &tmp_random_i,
                         int &tmp_random_j, 
                         arma::mat &data, 
                         int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size();
  // If the edge is already there, we need to substract one of the degrees
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  
  res.at(0) = -(degree_i == 0) -(degree_j == 0);
  res.at(1) = 0 ;
  res.at(2) = res.at(1) -res.at(0);
  return(res);
}

arma::vec stat_isolates_neg( const std::unordered_map< int, std::set<int>> &edges_pos,
                             const std::unordered_map< int, std::set<int>> &edges_neg,
                             int &n_actors,
                             int &tmp_random_i,
                             int &tmp_random_j, 
                             arma::mat &data, 
                             int &type){
  arma::vec res(3);
  int degree_i = edges_neg.at(tmp_random_i).size();
  int degree_j = edges_neg.at(tmp_random_j).size();
  // If the edge is already there, we need to substract one of the degrees
  if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  
  res.at(0) = 0;
  res.at(1) = -(degree_i == 0) -(degree_j == 0) ;
  res.at(2) = res.at(1) -res.at(0);
  return(res);
}
//  Calculate and return  
//  delta(present_val)^(0->1) 
//  delta(present_val)^(0->-1) 
//  delta(present_val)^(1->-1) 
// Dyadic Covariate neutral ->
arma::vec stat_cov_dyad_neut(const std::unordered_map< int, std::set<int>> &edges_pos,
                              const std::unordered_map< int, std::set<int>> &edges_neg,
                              int &n_actors,
                              int &tmp_random_i,
                              int &tmp_random_j, 
                              arma::mat &data, 
                              int &type){
  arma::vec res(3);
  res = {-data(tmp_random_i-1, tmp_random_j-1),-data(tmp_random_i-1, tmp_random_j-1),0};
  return(res);
}

// Dyadic Covariate negative
arma::vec stat_cov_dyad_pos(  const std::unordered_map< int, std::set<int>> &edges_pos,
                              const std::unordered_map< int, std::set<int>> &edges_neg,
                              int &n_actors,
                              int &tmp_random_i,
                              int &tmp_random_j, 
                              arma::mat &data, 
                              int &type){
  arma::vec res(3);
  res = {data(tmp_random_i-1, tmp_random_j-1),0,-data(tmp_random_i-1, tmp_random_j-1)};
  return(res);
}
// Dyadic Covariate negative
arma::vec stat_cov_dyad_neg(  const std::unordered_map< int, std::set<int>> &edges_pos,
                              const std::unordered_map< int, std::set<int>> &edges_neg,
                              int &n_actors,
                              int &tmp_random_i,
                              int &tmp_random_j, 
                              arma::mat &data, 
                              int &type){
  // Rcout << dyadic_data << std::endl;
  
  arma::vec res(3);
  res = {0,data(tmp_random_i-1, tmp_random_j-1),data(tmp_random_i-1, tmp_random_j-1)};
  return(res);
}

// Dyadic Covariate on both (negative and positiveS)
arma::vec stat_cov_dyad(  const std::unordered_map< int, std::set<int>> &edges_pos,
                          const std::unordered_map< int, std::set<int>> &edges_neg,
                          int &n_actors,
                          int &tmp_random_i,
                          int &tmp_random_j, 
                          arma::mat &data, 
                          int &type){
  // Rcout << dyadic_data << std::endl;
  
  arma::vec res(3);
  res = {data(tmp_random_i-1, tmp_random_j-1),data(tmp_random_i-1, tmp_random_j-1),0};
  return(res);
}

arma::vec stat_two_star_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
                            const std::unordered_map< int, std::set<int>> &edges_neg,
                            int &n_actors,
                            int &tmp_random_i,
                            int &tmp_random_j, 
                            arma::mat &data, 
                            int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size();
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    res.at(0) = degree_i + degree_j-2;
    res.at(1) = 0;
    res.at(2) = -res.at(0);
  } else {
    res.at(0) = degree_i + degree_j;
    res.at(1) = 0;
    res.at(2) = -res.at(0);
  }
  return(res);
}

arma::vec stat_two_star_neg(const std::unordered_map< int, std::set<int>> &edges_pos,
                            const std::unordered_map< int, std::set<int>> &edges_neg,
                            int &n_actors,
                            int &tmp_random_i,
                            int &tmp_random_j, 
                            arma::mat &data, 
                            int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size();
  
  
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    res.at(0) = 0;
    res.at(1) = degree_i + degree_j-2;
    res.at(2) = res.at(1);
  } else {
    res.at(0) = 0;
    res.at(1) = degree_i + degree_j;
    res.at(2) = res.at(1);
  }
  return(res);
}


void print_set(std::set<int> tmp){
  std::set<int>::iterator itr;
  for (itr = tmp.begin(); itr != tmp.end(); itr++) {
    Rcout << *itr << std::endl;
  }
}

arma::vec stat_gwdsp_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
                         const  std::unordered_map< int, std::set<int>> &edges_neg,
                         int &n_actors,
                         int &tmp_random_i,
                         int &tmp_random_j, 
                         arma::mat &data, 
                         int &type){
  arma::vec res(3);
  res.fill(0);
  double expo_min = (1-exp(-data.at(0,0)));  
  // double expo_pos = exp(data.at(0,0));  
  
  std::set<int> intersection, tmp,intersection_ij;
  // Step 1: Get all connections of tmp_random_i
  std::set<int> connections_of_i = edges_pos.at(tmp_random_i);
  // Step 2: Get all connections of tmp_random_j
  std::set<int> connections_of_j = edges_pos.at(tmp_random_j);
  // If i and j are already connected exclude them from the respective sets 
  connections_of_j.erase(tmp_random_i);
  connections_of_i.erase(tmp_random_j);
  // connections_of_j.insert(tmp_random_j);
  // connections_of_i.insert(tmp_random_i);
  // Step 3: Go through all connections of j and count its two-paths with i (call them h)
  std::set<int>::iterator itr;
  
  // Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  // Rcout << "Connections of i" << std::endl;
  // print_set(connections_of_i);
  // Rcout << "Connections of j" << std::endl;
  // print_set(connections_of_j);
  
  for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
    // For each connection of j we go through the positive edges 
    // and count how many of them are intersecting with connections of i 
    tmp = edges_pos.at(*itr);
    // Get the size of the intersection of tmp and connections_of_i
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_i), std::end(connections_of_i),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(0) +=  pow(expo_min, intersection.size());
    // Rcout << "Checking for "  + std::to_string(*itr) << std::endl;
    // Rcout << "Connections of tmp" << std::endl;
    // print_set(tmp);
    // Rcout << "Size" << std::endl;
    // Rcout <<intersection.size() << std::endl;
    
    
  }
  // Step 4: Go through all connections of i and count its two-paths with j
  for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
    // For each connection of j we go through the positive edges 
    tmp = edges_pos.at(*itr);
    // Get the size of the intersection of tmp and connections_of_i
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_j), std::end(connections_of_j),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(0) +=  pow(expo_min, intersection.size());
    // Rcout << "Connections of tmp" << std::endl;
    // print_set(tmp);
    // Rcout << "Size" << std::endl;
    // Rcout <<intersection.size() << std::endl;
  }
  // What to do with the intersection of connections of i and j 
  // std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
  //                       std::begin(connections_of_j), std::end(connections_of_j),
  //                       std::inserter(intersection_ij, std::begin(intersection_ij)));
  // if(intersection_ij.size()>0){
  //   // For each entry in there search with shared partners with i and j
  //   for (itr = intersection_ij.begin(); itr != intersection_ij.end(); itr++) {
  //     intersection.clear();
  //     std::set_intersection(std::begin(edges_pos.at(*itr)), std::end(edges_pos.at(*itr)),
  //                           std::begin(connections_of_j), std::end(connections_of_j),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     if(intersection.size()>0){
  //       res.at(0) +=  pow(expo_min, intersection.size()-1);
  //     }
  //     intersection.clear();
  //     std::set_intersection(std::begin(edges_pos.at(*itr)), std::end(edges_pos.at(*itr)),
  //                           std::begin(connections_of_i), std::end(connections_of_i),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     if(intersection.size()>0){
  //       res.at(0) +=  pow(expo_min, intersection.size()-1);
  //     }
  //   }
  // }
  res.at(1) = 0;
  res.at(2) = -res.at(0);
  
  // if(arma::sign(res.at(0)) == -1){
  //   Rcout << "Start" << std::endl;
  //   Rcout << res << std::endl;
  //   Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  //   Rcout << "Connections of i" << std::endl;
  //   print_set(connections_of_i);
  //   Rcout << "Connections of j" << std::endl;
  //   print_set(connections_of_j);
  //   
  //   arma::vec res(3);
  //   res.fill(0);
  //   double expo_min = (1-exp(-data.at(0,0)));  
  //   double expo_pos = exp(data.at(0,0));  
  //   Rcout << expo_min << std::endl;
  //   Rcout << data.at(0,0) << std::endl;
  //   
  //   
  //   std::set<int> intersection, tmp,intersection_ij;
  //   // Step 1: Get all connections of tmp_random_i
  //   std::set<int> connections_of_i = edges_pos.at(tmp_random_i);
  //   // Step 2: Get all connections of tmp_random_j
  //   std::set<int> connections_of_j = edges_pos.at(tmp_random_j);
  //   // If i and j are already connected exclude them from the respective sets 
  //   connections_of_j.erase(tmp_random_i);
  //   connections_of_i.erase(tmp_random_j);
  //   // connections_of_j.insert(tmp_random_j);
  //   // connections_of_i.insert(tmp_random_i);
  //   // Step 3: Go through all connections of j and count its two-paths with i (call them h)
  //   std::set<int>::iterator itr;
  //   
  //   // Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  //   // Rcout << "Connections of i" << std::endl;
  //   // print_set(connections_of_i);
  //   // Rcout << "Connections of j" << std::endl;
  //   // print_set(connections_of_j);
  //   
  //   for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
  //     // For each connection of j we go through the positive edges 
  //     // and count how many of them are intersecting with connections of i 
  //     tmp = edges_pos.at(*itr);
  //     // Get the size of the intersection of tmp and connections_of_i
  //     intersection.clear();    
  //     std::set_intersection(std::begin(tmp), std::end(tmp),
  //                           std::begin(connections_of_i), std::end(connections_of_i),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     res.at(0) +=  pow(expo_min, intersection.size());
  //     Rcout << "Result is "  + std::to_string(res.at(0)) << std::endl;
  //     
  //     Rcout << "Checking for "  + std::to_string(*itr) << std::endl;
  //     Rcout << "Connections of tmp" << std::endl;
  //     print_set(tmp);
  //     Rcout << "Size" << std::endl;
  //     Rcout <<intersection.size() << std::endl;
  //     
  //     
  //   }
  //   // Step 4: Go through all connections of i and count its two-paths with j
  //   for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
  //     // For each connection of j we go through the positive edges 
  //     tmp = edges_pos.at(*itr);
  //     // Get the size of the intersection of tmp and connections_of_i
  //     intersection.clear();    
  //     std::set_intersection(std::begin(tmp), std::end(tmp),
  //                           std::begin(connections_of_j), std::end(connections_of_j),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     res.at(0) +=  pow(expo_min, intersection.size());
  //     Rcout << "Result is "  + std::to_string(res.at(0)) << std::endl;
  //     
  //     Rcout << "Checking for "  + std::to_string(*itr) << std::endl;
  //     Rcout << "Connections of tmp" << std::endl;
  //     print_set(tmp);
  //     Rcout << "Size" << std::endl;
  //     Rcout <<intersection.size() << std::endl;
  //   }
  // }
  
  return(res);
}

arma::vec stat_gwdsp_neg(const std::unordered_map< int, std::set<int>> &edges_pos,
                         const  std::unordered_map< int, std::set<int>> &edges_neg,
                         int &n_actors,
                         int &tmp_random_i,
                         int &tmp_random_j, 
                         arma::mat &data, 
                         int &type){
  arma::vec res(3);
  res.fill(0);
  double expo_min = (1-exp(-data.at(0,0)));  
  // double expo_pos = exp(data.at(0,0));  
  
  std::set<int> intersection, tmp,intersection_ij;
  // Step 1: Get all connections of tmp_random_i
  std::set<int> connections_of_i = edges_neg.at(tmp_random_i);
  // Step 2: Get all connections of tmp_random_j
  std::set<int> connections_of_j = edges_neg.at(tmp_random_j);
  // If i and j are already connected exclude them from the respective sets 
  connections_of_j.erase(tmp_random_i);
  connections_of_i.erase(tmp_random_j);
  // connections_of_j.insert(tmp_random_j);
  // connections_of_i.insert(tmp_random_i);
  // Step 3: Go through all connections of j and count its two-paths with i (call them h)
  std::set<int>::iterator itr;
  
  // Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  // Rcout << "Connections of i" << std::endl;
  // print_set(connections_of_i);
  // Rcout << "Connections of j" << std::endl;
  // print_set(connections_of_j);
  
  for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
    // For each connection of j we go through the positive edges 
    // and count how many of them are intersecting with connections of i 
    tmp = edges_neg.at(*itr);
    // Get the size of the intersection of tmp and connections_of_i
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_i), std::end(connections_of_i),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(1) +=  pow(expo_min, intersection.size());
    // Rcout << "Checking for "  + std::to_string(*itr) << std::endl;
    // Rcout << "Connections of tmp" << std::endl;
    // print_set(tmp);
    // Rcout << "Size" << std::endl;
    // Rcout <<intersection.size() << std::endl;
    
    
  }
  // Step 4: Go through all connections of i and count its two-paths with j
  for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
    // For each connection of j we go through the positive edges 
    tmp = edges_neg.at(*itr);
    // Get the size of the intersection of tmp and connections_of_i
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_j), std::end(connections_of_j),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(1) +=  pow(expo_min, intersection.size());
    // Rcout << "Connections of tmp" << std::endl;
    // print_set(tmp);
    // Rcout << "Size" << std::endl;
    // Rcout <<intersection.size() << std::endl;
  }
  // What to do with the intersection of connections of i and j 
  // std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
  //                       std::begin(connections_of_j), std::end(connections_of_j),
  //                       std::inserter(intersection_ij, std::begin(intersection_ij)));
  // if(intersection_ij.size()>0){
  //   // For each entry in there search with shared partners with i and j
  //   for (itr = intersection_ij.begin(); itr != intersection_ij.end(); itr++) {
  //     intersection.clear();
  //     std::set_intersection(std::begin(edges_pos.at(*itr)), std::end(edges_pos.at(*itr)),
  //                           std::begin(connections_of_j), std::end(connections_of_j),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     if(intersection.size()>0){
  //       res.at(0) +=  pow(expo_min, intersection.size()-1);
  //     }
  //     intersection.clear();
  //     std::set_intersection(std::begin(edges_pos.at(*itr)), std::end(edges_pos.at(*itr)),
  //                           std::begin(connections_of_i), std::end(connections_of_i),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     if(intersection.size()>0){
  //       res.at(0) +=  pow(expo_min, intersection.size()-1);
  //     }
  //   }
  // }
  res.at(0) = 0;
  res.at(2) = res.at(1);
  
  // if(arma::sign(res.at(0)) == -1){
  //   Rcout << "Start" << std::endl;
  //   Rcout << res << std::endl;
  //   Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  //   Rcout << "Connections of i" << std::endl;
  //   print_set(connections_of_i);
  //   Rcout << "Connections of j" << std::endl;
  //   print_set(connections_of_j);
  //   
  //   arma::vec res(3);
  //   res.fill(0);
  //   double expo_min = (1-exp(-data.at(0,0)));  
  //   double expo_pos = exp(data.at(0,0));  
  //   Rcout << expo_min << std::endl;
  //   Rcout << data.at(0,0) << std::endl;
  //   
  //   
  //   std::set<int> intersection, tmp,intersection_ij;
  //   // Step 1: Get all connections of tmp_random_i
  //   std::set<int> connections_of_i = edges_pos.at(tmp_random_i);
  //   // Step 2: Get all connections of tmp_random_j
  //   std::set<int> connections_of_j = edges_pos.at(tmp_random_j);
  //   // If i and j are already connected exclude them from the respective sets 
  //   connections_of_j.erase(tmp_random_i);
  //   connections_of_i.erase(tmp_random_j);
  //   // connections_of_j.insert(tmp_random_j);
  //   // connections_of_i.insert(tmp_random_i);
  //   // Step 3: Go through all connections of j and count its two-paths with i (call them h)
  //   std::set<int>::iterator itr;
  //   
  //   // Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  //   // Rcout << "Connections of i" << std::endl;
  //   // print_set(connections_of_i);
  //   // Rcout << "Connections of j" << std::endl;
  //   // print_set(connections_of_j);
  //   
  //   for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
  //     // For each connection of j we go through the positive edges 
  //     // and count how many of them are intersecting with connections of i 
  //     tmp = edges_pos.at(*itr);
  //     // Get the size of the intersection of tmp and connections_of_i
  //     intersection.clear();    
  //     std::set_intersection(std::begin(tmp), std::end(tmp),
  //                           std::begin(connections_of_i), std::end(connections_of_i),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     res.at(0) +=  pow(expo_min, intersection.size());
  //     Rcout << "Result is "  + std::to_string(res.at(0)) << std::endl;
  //     
  //     Rcout << "Checking for "  + std::to_string(*itr) << std::endl;
  //     Rcout << "Connections of tmp" << std::endl;
  //     print_set(tmp);
  //     Rcout << "Size" << std::endl;
  //     Rcout <<intersection.size() << std::endl;
  //     
  //     
  //   }
  //   // Step 4: Go through all connections of i and count its two-paths with j
  //   for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
  //     // For each connection of j we go through the positive edges 
  //     tmp = edges_pos.at(*itr);
  //     // Get the size of the intersection of tmp and connections_of_i
  //     intersection.clear();    
  //     std::set_intersection(std::begin(tmp), std::end(tmp),
  //                           std::begin(connections_of_j), std::end(connections_of_j),
  //                           std::inserter(intersection, std::begin(intersection)));
  //     res.at(0) +=  pow(expo_min, intersection.size());
  //     Rcout << "Result is "  + std::to_string(res.at(0)) << std::endl;
  //     
  //     Rcout << "Checking for "  + std::to_string(*itr) << std::endl;
  //     Rcout << "Connections of tmp" << std::endl;
  //     print_set(tmp);
  //     Rcout << "Size" << std::endl;
  //     Rcout <<intersection.size() << std::endl;
  //   }
  // }
  
  return(res);
}

arma::vec stat_gwdasp(const std::unordered_map< int, std::set<int>> &edges_pos,
                      const  std::unordered_map< int, std::set<int>> &edges_neg,
                      int &n_actors,
                      int &tmp_random_i,
                      int &tmp_random_j, 
                      arma::mat &data, 
                      int &type){
  arma::vec res(3);
  res.fill(0);
  double expo_min = (1-exp(-data.at(0,0)));  
  // double expo_pos = exp(data.at(0,0));  
  
  std::set<int> intersection, tmp,intersection_ij;
  std::set<int> connections_of_i_neg = edges_neg.at(tmp_random_i);
  std::set<int> connections_of_j_neg = edges_neg.at(tmp_random_j);
  std::set<int> connections_of_i_pos = edges_pos.at(tmp_random_i);
  std::set<int> connections_of_j_pos = edges_pos.at(tmp_random_j);
  // If i and j are already connected exclude them from the respective sets 
  connections_of_j_pos.erase(tmp_random_i);
  connections_of_i_pos.erase(tmp_random_j);
  connections_of_j_neg.erase(tmp_random_i);
  connections_of_i_neg.erase(tmp_random_j);
  // Step 3: Go through all connections of j and count its two-paths with i (call them h)
  std::set<int>::iterator itr;
  
  // Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  // Rcout << "Connections of i" << std::endl;
  // print_set(connections_of_i);
  // Rcout << "Connections of j" << std::endl;
  // print_set(connections_of_j);
  
  for (itr = connections_of_j_neg.begin(); itr != connections_of_j_neg.end(); itr++) {
    tmp = edges_neg.at(*itr);
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_i_pos), std::end(connections_of_i_pos),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(0) +=  pow(expo_min, intersection.size());
  }
  for (itr = connections_of_i_neg.begin(); itr != connections_of_i_neg.end(); itr++) {
    tmp = edges_neg.at(*itr);
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_j_pos), std::end(connections_of_j_pos),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(0) +=  pow(expo_min, intersection.size());
  }
  // Do it for 0 to -1
  for (itr = connections_of_j_pos.begin(); itr != connections_of_j_pos.end(); itr++) {
    tmp = edges_pos.at(*itr);
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_i_neg), std::end(connections_of_i_neg),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(1) +=  pow(expo_min, intersection.size());
  }
  for (itr = connections_of_i_pos.begin(); itr != connections_of_i_pos.end(); itr++) {
    tmp = edges_pos.at(*itr);
    intersection.clear();    
    std::set_intersection(std::begin(tmp), std::end(tmp),
                          std::begin(connections_of_j_neg), std::end(connections_of_j_neg),
                          std::inserter(intersection, std::begin(intersection)));
    res.at(1) +=  pow(expo_min, intersection.size());
  }
  res.at(2) = -res.at(0) + res.at(1);
  
  return(res);
}


arma::vec stat_gwesf_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
                         const std::unordered_map< int, std::set<int>> &edges_neg,
                         int &n_actors,
                         int &tmp_random_i,
                         int &tmp_random_j, 
                         arma::mat &data, 
                         int &type){
  arma::vec res(3);
  double expo_min = (1-exp(-data.at(0,0)));  
  double expo_pos = exp(data.at(0,0));  
  
  std::set<int>::iterator itr;
  std::set<int> intersection_ij,intersection, tmp;
  // Step 1: Get all positive connections of tmp_random_i
  std::set<int> connections_of_i = edges_pos.at(tmp_random_i);
  // Step 2: Get all positive connections of tmp_random_j
  std::set<int> connections_of_j = edges_pos.at(tmp_random_j);
  // If i and j are already connected exclude them from the respective sets 
  connections_of_j.erase(tmp_random_i);
  connections_of_i.erase(tmp_random_j);
  // Step 3: Get the intersection of connections to i and j -> case a)
  std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
                        std::begin(connections_of_j), std::end(connections_of_j),
                        std::inserter(intersection_ij, std::begin(intersection_ij)));
  if(intersection_ij.size()>0){
    res.at(0) += expo_pos*(1- pow(expo_min, intersection_ij.size()));
    connections_of_j.insert(tmp_random_j);
    connections_of_i.insert(tmp_random_i);
    // Step 4.1: Get all connections from actors in this intersection 
    for (itr = intersection_ij.begin(); itr != intersection_ij.end(); itr++) {
      tmp = edges_pos.at(*itr);
      // Rcout << "Connections of tmp" << std::endl;
      // print_set(tmp);
      // Rcout << "Connections of i" << std::endl;
      // print_set(connections_of_i);
      // Rcout << "Connections of j" << std::endl;
      // print_set(connections_of_j);
      // Get the size of the intersection of tmp and connections_of_i -> case b)
      intersection.clear();    
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_i), std::end(connections_of_i),
                            std::inserter(intersection, std::begin(intersection)));
      if(intersection.size()!=0){
        res.at(0) +=  pow(expo_min, intersection.size() -1);
      }
      // Get the size of the intersection of tmp and connections_of_j -> case c)
      // Rcout << intersection.size() << std::endl;
      
      intersection.clear();  
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_j), std::end(connections_of_j),
                            std::inserter(intersection, std::begin(intersection)));
      if(intersection.size()!=0){
        res.at(0) +=  pow(expo_min, intersection.size()-1);
      }
      // Rcout << intersection.size() << std::endl;
      tmp.clear();
    }
  }
  res.at(1) = 0;
  res.at(2) = -res.at(0);
  // Rcout << res << std::endl;
  return(res);
}

arma::vec stat_gwese_neg(const std::unordered_map< int, std::set<int>> &edges_pos,
                         const std::unordered_map< int, std::set<int>> &edges_neg,
                         int &n_actors,
                         int &tmp_random_i,
                         int &tmp_random_j, 
                         arma::mat &data, 
                         int &type){
  arma::vec res(3);
  double expo_min = (1-exp(-data.at(0,0)));  
  double expo_pos = exp(data.at(0,0));  
  
  std::set<int>::iterator itr;
  std::set<int> intersection_ij,intersection, tmp;
  // Step 1: Get all positive connections of tmp_random_i
  std::set<int> connections_of_i = edges_neg.at(tmp_random_i);
  // Step 2: Get all positive connections of tmp_random_j
  std::set<int> connections_of_j = edges_neg.at(tmp_random_j);
  // If i and j are already connected exclude them from the respective sets 
  connections_of_j.erase(tmp_random_i);
  connections_of_i.erase(tmp_random_j);
  // Step 3: Get the intersection of connections to i and j -> case a)
  std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
                        std::begin(connections_of_j), std::end(connections_of_j),
                        std::inserter(intersection_ij, std::begin(intersection_ij)));
  if(intersection_ij.size()>0){
    res.at(1) += expo_pos*(1- pow(expo_min, intersection_ij.size()));
    connections_of_j.insert(tmp_random_j);
    connections_of_i.insert(tmp_random_i);
    // Step 4.1: Get all connections from actors in this intersection 
    for (itr = intersection_ij.begin(); itr != intersection_ij.end(); itr++) {
      tmp = edges_neg.at(*itr);
      // Rcout << "Connections of tmp" << std::endl;
      // print_set(tmp);
      // Rcout << "Connections of i" << std::endl;
      // print_set(connections_of_i);
      // Rcout << "Connections of j" << std::endl;
      // print_set(connections_of_j);
      // Get the size of the intersection of tmp and connections_of_i -> case b)
      intersection.clear();    
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_i), std::end(connections_of_i),
                            std::inserter(intersection, std::begin(intersection)));
      if(intersection.size()!=0){
        res.at(1) +=  pow(expo_min, intersection.size() -1);
      }
      // Get the size of the intersection of tmp and connections_of_j -> case c)
      // Rcout << intersection.size() << std::endl;
      
      intersection.clear();  
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_j), std::end(connections_of_j),
                            std::inserter(intersection, std::begin(intersection)));
      if(intersection.size()!=0){
        res.at(1) +=  pow(expo_min, intersection.size()-1);
      }
      // Rcout << intersection.size() << std::endl;
      tmp.clear();
    }
  }
  res.at(0) = 0;
  res.at(2) = res.at(1);
  // Rcout << res << std::endl;
  return(res);
}

arma::vec stat_gwese_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
                     const std::unordered_map< int, std::set<int>> &edges_neg,
                     int &n_actors,
                     int &tmp_random_i,
                     int &tmp_random_j, 
                     arma::mat &data, 
                     int &type){
  arma::vec res(3);
  res.fill(0);
  double expo_min = (1-exp(-data.at(0,0)));  
  double expo_pos = exp(data.at(0,0));  
  
  std::set<int>::iterator itr;
  std::set<int> intersection_ij,intersection, tmp;
  std::set<int> connections_of_i_neg = edges_neg.at(tmp_random_i);
  std::set<int> connections_of_j_neg = edges_neg.at(tmp_random_j);
  std::set<int> connections_of_i_pos = edges_pos.at(tmp_random_i);
  std::set<int> connections_of_j_pos = edges_pos.at(tmp_random_j);
  // Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  // 
  // Rcout << "Positive connections of i" << std::endl;
  // print_set(connections_of_i_pos);
  // Rcout << "Positive connections of j" << std::endl;
  // print_set(connections_of_j_pos);
  // Rcout << "Negative connections of i" << std::endl;
  // print_set(connections_of_i_neg);
  // Rcout << "Negative connections of j" << std::endl;
  // print_set(connections_of_j_neg);
  // 
  // Step 3: Get the intersection of connections to i and j -> case a)
  // Rcout << "Positive connections of i" << std::endl;
  std::set_intersection(std::begin(connections_of_i_neg), std::end(connections_of_i_neg),
                        std::begin(connections_of_j_neg), std::end(connections_of_j_neg),
                        std::inserter(intersection_ij, std::begin(intersection_ij)));
  res.at(0) += expo_pos*(1- pow(expo_min, intersection_ij.size()));
  
  // If tmp_i and tmp_j are already connected, we need to delete them from the tmp graph 
  connections_of_i_neg.erase(tmp_random_j);
  connections_of_j_neg.erase(tmp_random_i);
  // Rcout << "Positive connections of i" << std::endl;
  // Rcout <<intersection_ij.size() << std::endl;
  // connections_of_j_neg.insert(tmp_random_j);
  // connections_of_i_neg.insert(tmp_random_i);
  // Rcout << "Positive connections of i" << std::endl;
  // if(connections_of_i_pos.size()>0){
  //   
  // }
  // Check for each positive connection of i that has a negative connection to j how many negative common partner it has with i
  // How many common enemies do we have?  
  for (itr = connections_of_i_pos.begin(); itr != connections_of_i_pos.end(); itr++) {
    tmp = edges_neg.at(*itr);
    if(tmp.count(tmp_random_j)){
      intersection.clear();    
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_i_neg), std::end(connections_of_i_neg),
                            std::inserter(intersection, std::begin(intersection)));
      res.at(1) +=  pow(expo_min, intersection.size());
      // Rcout << intersection.size() << std::endl;
      
    }
  }
  // Rcout << "Positive connections of i" << std::endl;
  // Check for each positive connection of j that has a negative connection to i how many negative common partner it has with j
  // How many common enemies do we have?  
  for (itr = connections_of_j_pos.begin(); itr != connections_of_j_pos.end(); itr++) {
    tmp = edges_neg.at(*itr);
    if(tmp.count(tmp_random_i)){
      intersection.clear();    
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_j_neg), std::end(connections_of_j_neg),
                            std::inserter(intersection, std::begin(intersection)));
      res.at(1) +=  pow(expo_min, intersection.size());
      
      // Rcout << intersection.size() << std::endl;
    }
  }
  // Rcout << "Positive connections of i" << std::endl;
  res.at(2) = res.at(1) - res.at(0);
  // Rcout << res << std::endl;
  return(res);
}

arma::vec stat_gwesf_neg(const std::unordered_map< int, std::set<int>> &edges_pos,
                     const std::unordered_map< int, std::set<int>> &edges_neg,
                     int &n_actors,
                     int &tmp_random_i,
                     int &tmp_random_j, 
                     arma::mat &data, 
                     int &type){
  arma::vec res(3);
  res.fill(0);
  double expo_min = (1-exp(-data.at(0,0)));  
  double expo_pos = exp(data.at(0,0));  
  
  std::set<int>::iterator itr;
  std::set<int> intersection_ij,intersection, tmp;
  std::set<int> connections_of_i_neg = edges_neg.at(tmp_random_i);
  std::set<int> connections_of_j_neg = edges_neg.at(tmp_random_j);
  std::set<int> connections_of_i_pos = edges_pos.at(tmp_random_i);
  std::set<int> connections_of_j_pos = edges_pos.at(tmp_random_j);
  // Rcout << "Between " + std::to_string(tmp_random_i) +  " and "  + std::to_string(tmp_random_j) << std::endl;
  // 
  // Rcout << "Positive connections of i" << std::endl;
  // print_set(connections_of_i_pos);
  // Rcout << "Positive connections of j" << std::endl;
  // print_set(connections_of_j_pos);
  // Rcout << "Negative connections of i" << std::endl;
  // print_set(connections_of_i_neg);
  // Rcout << "Negative connections of j" << std::endl;
  // print_set(connections_of_j_neg);
  // 
  // Step 3: Get the intersection of connections to i and j -> case a)
  // Rcout << "Positive connections of i" << std::endl;
  std::set_intersection(std::begin(connections_of_i_pos), std::end(connections_of_i_pos),
                        std::begin(connections_of_j_pos), std::end(connections_of_j_pos),
                        std::inserter(intersection_ij, std::begin(intersection_ij)));
  res.at(1) += expo_pos*(1- pow(expo_min, intersection_ij.size()));
  
  // If tmp_i and tmp_j are already connected, we need to delete them from the tmp graph 
  connections_of_i_pos.erase(tmp_random_j);
  connections_of_j_pos.erase(tmp_random_i);
  // Rcout << "Positive connections of i" << std::endl;
  // Rcout <<intersection_ij.size() << std::endl;
  // connections_of_j_neg.insert(tmp_random_j);
  // connections_of_i_neg.insert(tmp_random_i);
  // Rcout << "Positive connections of i" << std::endl;
  // if(connections_of_i_pos.size()>0){
  //   
  // }
  // Check for each negative connection of i that has a negative connection to j how many negative common partner it has with i
  // How many common friends do we have?  
  for (itr = connections_of_i_neg.begin(); itr != connections_of_i_neg.end(); itr++) {
    tmp = edges_pos.at(*itr);
    if(tmp.count(tmp_random_j)){
      intersection.clear();    
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_i_pos), std::end(connections_of_i_pos),
                            std::inserter(intersection, std::begin(intersection)));
      res.at(0) +=  pow(expo_min, intersection.size());
      // Rcout << intersection.size() << std::endl;
      
    }
  }
  // Rcout << "Positive connections of i" << std::endl;
  // Check for each positive connection of j that has a negative connection to i how many negative common partner it has with j
  // How many common friends do we have?  
  for (itr = connections_of_j_neg.begin(); itr != connections_of_j_neg.end(); itr++) {
    tmp = edges_pos.at(*itr);
    if(tmp.count(tmp_random_i)){
      intersection.clear();    
      std::set_intersection(std::begin(tmp), std::end(tmp),
                            std::begin(connections_of_j_pos), std::end(connections_of_j_pos),
                            std::inserter(intersection, std::begin(intersection)));
      res.at(0) +=  pow(expo_min, intersection.size());
      
      // Rcout << intersection.size() << std::endl;
    } 
  }
  // Rcout << "Positive connections of i" << std::endl;
  res.at(2) = res.at(1) - res.at(0);
  // Rcout << res << std::endl;
  return(res);
}

arma::vec stat_gwesp_pos_old(const std::unordered_map< int, std::set<int>> &edges_pos,
                             const std::unordered_map< int, std::set<int>> &edges_neg,
                             int &n_actors,
                             int &tmp_random_i,
                             int &tmp_random_j, 
                             arma::mat &data, 
                             int &type){
  arma::vec res(3);
  std::set<int> tmp;
  double expo_min = (1-exp(-data.at(0,0)));  
  double expo_pos = exp(data.at(0,0));  
  // Step 1: Get all connections of tmp_random_i
  std::set<int> connections_of_i = edges_pos.at(tmp_random_i);
  // Step 2: Get all connections of tmp_random_j
  std::set<int> connections_of_j = edges_pos.at(tmp_random_j);
  // If i and j are already connected exclude them from the respective sets 
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    connections_of_j.erase(tmp_random_i);
    connections_of_i.erase(tmp_random_j);
  } 
  // Step 3: Get all connections of connections of tmp_random_j
  std::vector<std::set<int>> edges_to_j;
  std::vector<std::set<int>> edges_to_i;
  
  std::set<int>::iterator itr;
  for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
    // Rcout << "Accessing " + std::to_string(*itr) << std::endl;
    edges_to_j.push_back(edges_pos.at(*itr));
  } 
  for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
    // Rcout << "Accessing " + std::to_string(*itr) << std::endl;
    edges_to_i.push_back(edges_pos.at(*itr));
    // edges_to_i[n] = edges_pos.at(*itr);
  }
  
  // Step 4: 
  // Rcout << "Between " + std::to_string(tmp_random_i) + " and " + std::to_string(tmp_random_j)<< std::endl;
  // Rcout <<  connections_of_i.size() << std::endl;
  // Rcout <<  edges_to_j.size() << std::endl;
  
  for(unsigned int i = 0; i<(edges_to_j.size()); i +=1){
    // Rcout <<  "Look at " + std::to_string(i) << std::endl;
    tmp.clear();    
    std::set_intersection(std::begin(edges_to_j.at(i)), std::end(edges_to_j.at(i)),
                          std::begin(connections_of_i), std::end(connections_of_i),
                          std::inserter(tmp, std::begin(tmp)));
    if(tmp.size()>0){
      // res.at(0) += exp(-data.at(0,0)*(tmp.size())) - exp(-data.at(0,0)*(tmp.size()-1));
      res.at(0) +=  expo_pos*(pow(expo_min, tmp.size()) - pow(expo_min, tmp.size()+1));
    }
  }
  
  for(unsigned int i = 0; i<(edges_to_i.size()); i +=1){
    // Rcout <<  "Look at " + std::to_string(i) << std::endl;
    // Rcout <<  "From size " + std::to_string(edges_to_j.size()) << std::endl;
    tmp.clear();    
    
    std::set_intersection(std::begin(edges_to_i.at(i)), std::end(edges_to_i.at(i)),
                          std::begin(connections_of_j), std::end(connections_of_j),
                          std::inserter(tmp, std::begin(tmp)));
    if(tmp.size()>0){
      res.at(0) +=  expo_pos*(pow(expo_min, tmp.size()) - pow(expo_min, tmp.size()+1));
    } 
  } 
  res.at(1) = 0;
  res.at(2) = -res.at(0);
  return(res);
} 

// 
// arma::vec stat_gwesp_old(std::unordered_map< int, std::set<int>> &edges_pos,
//                          std::unordered_map< int, std::set<int>> &edges_neg,
//                          int &n_actors,
//                          int &tmp_random_i,
//                          int &tmp_random_j, 
//                          arma::mat &data, 
//                          int &type){
//   arma::vec res(3);
//   std::set<int> tmp;
//   // Step 1: Get all connections of tmp_random_i
//   std::set<int> connections_of_i = edges_pos.at(tmp_random_i);
//   // Step 2: Get all connections of tmp_random_j
//   std::set<int> connections_of_j = edges_pos.at(tmp_random_j);
//   // Step 3: Get all connections of connections of tmp_random_j
//   std::unordered_map< int, std::set<int>> edges_to_j;
//   std::unordered_map< int, std::set<int>> edges_to_i;
//   
//   std::set<int>::iterator itr;
//   int n = 1;
//   // Rcout << "Start" << std::endl;
//   
//   if(connections_of_j.size()>0){
//     for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
//       // Rcout << "Accessing " + std::to_string(*itr) << std::endl;
//       edges_to_j[n] = edges_pos.at(*itr);
//       n +=1;
//     }
//   }
//   n = 1;
//   if(connections_of_i.size()>0){
//     for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
//       // Rcout << "Accessing " + std::to_string(*itr) << std::endl;
//       edges_to_i[n] = edges_pos.at(*itr);
//       n +=1;
//     }
//   }
//   
//   // Step 4: 
//   // Rcout << "Between " + std::to_string(tmp_random_i) + " and " + std::to_string(tmp_random_j)<< std::endl;
//   // Rcout <<  connections_of_i.size() << std::endl;
//   // Rcout <<  edges_to_j.size() << std::endl;
//   
//   if(connections_of_i.size()>0 & edges_to_j.size()>0){
//     for(int i = 1; i<(edges_to_j.size()); i +=1){
//       // Rcout <<  "Look at " + std::to_string(i) << std::endl;
//       // Rcout <<  "From size " + std::to_string(edges_to_j.size()) << std::endl;
//       
//       std::set_intersection(std::begin(edges_to_j.at(i)), std::end(edges_to_j.at(i)),
//                             std::begin(connections_of_i), std::end(connections_of_i),
//                             std::inserter(tmp, std::begin(tmp)));
//       if(tmp.size()>0){
//         res.at(0) += exp(-data.at(0,0)*(tmp.size()-1)) - exp(-data.at(0,0)*tmp.size());
//       }
//     }
//   }
//   
//   if(connections_of_j.size()>0 & edges_to_i.size()>0){
//     for(int i = 1; i<(edges_to_i.size()); i +=1){
//       // Rcout <<  "Look at " + std::to_string(i) << std::endl;
//       // Rcout <<  "From size " + std::to_string(edges_to_j.size()) << std::endl;
//       
//       std::set_intersection(std::begin(edges_to_i.at(i)), std::end(edges_to_i.at(i)),
//                             std::begin(connections_of_j), std::end(connections_of_j),
//                             std::inserter(tmp, std::begin(tmp)));
//       if(tmp.size()>0){
//         res.at(0) += exp(-data.at(0,0)*(tmp.size()-1)) - exp(-data.at(0,0)*tmp.size());
//       }
//     }
//   }
//   res.at(1) = 0;
//   res.at(2) = -res.at(0);
//   return(res);
// }

arma::vec stat_gwdegree(const std::unordered_map< int, std::set<int>> &edges_pos,
                        const std::unordered_map< int, std::set<int>> &edges_neg,
                        int &n_actors,
                        int &tmp_random_i,
                        int &tmp_random_j, 
                        arma::mat &data, 
                        int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size() + edges_neg.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size() + edges_neg.at(tmp_random_j).size();
  double expo_min = (1-exp(-data.at(0,0)));  
  // double expo_pos = exp(data.at(0,0)); 
  // If the edge is already there, we need to substract one of the degrees
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  // Rcout << -(1-exp(-data.at(0,0)))*(exp(-data.at(0,0)*degree_i) + exp(-data.at(0,0)*degree_j)) << std::endl;
  
  res.at(0) = pow(expo_min, degree_i) + pow(expo_min, degree_j);
  res.at(1) = res.at(0);
  res.at(2) = 0;
  // res.at(2) = -(degree_i == type) -(degree_j == type) + (degree_i == (type+1))  + (degree_j == (type+1));
  return(res);
}

arma::vec stat_gwdegree_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
                            const std::unordered_map< int, std::set<int>> &edges_neg,
                            int &n_actors,
                            int &tmp_random_i,
                            int &tmp_random_j, 
                            arma::mat &data, 
                            int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size();
  double expo_min = (1-exp(-data.at(0,0)));  
  // double expo_pos = exp(data.at(0,0)); 
  // If the edge is already there, we need to substract one of the degrees
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  
  // Rcout << -(1-exp(-data.at(0,0)))*(exp(-data.at(0,0)*degree_i) + exp(-data.at(0,0)*degree_j)) << std::endl;
  
  res.at(0) = pow(expo_min, degree_i) + pow(expo_min, degree_j);
  res.at(1) = 0;
  res.at(2) = -res.at(0);
  // res.at(2) = -(degree_i == type) -(degree_j == type) + (degree_i == (type+1))  + (degree_j == (type+1));
  return(res);
}

arma::vec stat_gwdegree_neg(const std::unordered_map< int, std::set<int>> &edges_pos,
                            const std::unordered_map< int, std::set<int>> &edges_neg,
                            int &n_actors,
                            int &tmp_random_i,
                            int &tmp_random_j, 
                            arma::mat &data, 
                            int &type){
  arma::vec res(3);
  double expo_min = (1-exp(-data.at(0,0)));  
  // double expo_pos = exp(data.at(0,0));  
  int degree_i = edges_neg.at(tmp_random_i).size();
  int degree_j = edges_neg.at(tmp_random_j).size();
  // If the edge is already there, we need to substract one of the degrees
  if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  res.at(0) = 0;
  res.at(1) = pow(expo_min, degree_i) + pow(expo_min, degree_j);
  // res.at(1) = expo_pos*(pow(expo_min, degree_i) - pow(expo_min, degree_i+1)) + expo_pos*(pow(expo_min, degree_j) - pow(expo_min, degree_j+1));
  // res.at(1) = -(1-exp(-data.at(0,0)))*(exp(-data.at(0,0)*degree_i) + exp(-data.at(0,0)*degree_j));
  res.at(2) = res.at(1);
  return(res);
}

arma::vec stat_degree(const std::unordered_map< int, std::set<int>> &edges_pos,
                      const  std::unordered_map< int, std::set<int>> &edges_neg,
                      int &n_actors,
                      int &tmp_random_i,
                      int &tmp_random_j, 
                      arma::mat &data, 
                      int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size() + edges_neg.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size() + edges_neg.at(tmp_random_j).size();
  // If the edge is already there, we need to substract one of the degrees
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  // Rcout << degree_i << std::endl;
  // Rcout << degree_j << std::endl;
  
  // First entry
  // If deg_pos == type the number of actors with deg_pos == type will get smaller 
  // If deg_pos == type-1 the number of actors with deg_pos == type will get bigger 
  res.at(0) = -(degree_i == type) -(degree_j == type) + (degree_i == (type-1))  + (degree_j == (type-1));
  res.at(1) = res.at(0) ;
  res.at(2) = 0;
  return(res);
}

arma::vec stat_degree_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
                          const  std::unordered_map< int, std::set<int>> &edges_neg,
                          int &n_actors,
                          int &tmp_random_i,
                          int &tmp_random_j, 
                          arma::mat &data, 
                          int &type){
  arma::vec res(3);
  int degree_i = edges_pos.at(tmp_random_i).size();
  int degree_j = edges_pos.at(tmp_random_j).size();
  // If the edge is already there, we need to substract one of the degrees
  if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
    degree_i -= 1;
    degree_j -= 1;
  }
  // Rcout << degree_i << std::endl;
  // Rcout << degree_j << std::endl;
  
  // First entry
  // If deg_pos == type the number of actors with deg_pos == type will get smaller 
  // If deg_pos == type-1 the number of actors with deg_pos == type will get bigger 
  res.at(0) = -(degree_i == type) -(degree_j == type) + (degree_i == (type-1))  + (degree_j == (type-1));
  res.at(1) = 0;
  res.at(2) = -res.at(0);
  return(res);
}

arma::vec stat_degree_neg(const std::unordered_map< int, std::set<int>> &edges_pos,
                          const std::unordered_map< int, std::set<int>> &edges_neg,
                          int &n_actors,
                          int &tmp_random_i,
                          int &tmp_random_j, 
                          arma::mat &data, 
                          int &type){
  arma::vec res(3);
  int degree_i = edges_neg.at(tmp_random_i).size();
  int degree_j = edges_neg.at(tmp_random_j).size();
  // First entry
  // If deg_neg == type the number of actors with deg_neg == type will get smaller 
  // If deg_neg == type-1 the number of actors with deg_neg == type will get bigger 
  res.at(0) = 0;
  res.at(1) = -(degree_i == type) -(degree_j == type) + (degree_i == (type-1))  + (degree_j == (type-1));
  res.at(2) = res.at(1);
  return(res);
}

