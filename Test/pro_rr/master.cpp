#include "include/globals.h"
#include "include/setup.hpp"
#include "include/evolution.hpp"

// [[Rcpp::export]]
void roadrunner_setup(std::string antimony_file) {
    rr_setup(antimony_file);
}

// [[Rcpp::export]]
NumericVector roadrunner_evolution(NumericVector parent_state, Environment node) {
    return rr_evolution(parent_state, node);
}
