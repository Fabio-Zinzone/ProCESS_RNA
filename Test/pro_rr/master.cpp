#include "include/globals.h"
#include "include/setup.hpp"
#include "include/evolution.hpp"
#include "include/matrix.hpp"

using namespace Rcpp;

std::random_device rd;
std::mt19937 gen(rd());

// [[Rcpp::export]]
void set_cpp_seed(int seed) {
  gen.seed(seed);
}

// [[Rcpp::export]]
void roadrunner_setup(std::string antimony_file) {
  rr_setup(antimony_file);
}

// [[Rcpp::export]]
NumericVector roadrunner_evolution(NumericVector parent_state, Environment node) {
  return rr_evolution(parent_state, node);
}

// [[Rcpp::export]]
void apply_poisson() {
  int_poisson();
}

// [[Rcpp::export]]
Matrix_type save_grn_matrix() {
  return grn_matrix();
}

// [[Rcpp::export]]
void save_csv_matrix(std::string filename) {
  save_matrix(filename);
}
