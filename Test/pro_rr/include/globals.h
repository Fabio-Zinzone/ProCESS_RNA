#ifndef GLOBALS_H
#define GLOBALS_H

#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <fstream>
#include <Rcpp.h>
#include <rr/rrRoadRunner.h>
#include <rr/rrExecutableModel.h>
#include "antimony_api.h"

using namespace rr;
using Matrix_type = Rcpp::NumericMatrix;
//using Matrix_type = Rcpp::IntegerMatrix;

extern std::mt19937 gen;
extern RoadRunner* rrg;
extern std::map<int, std::vector<double>> grn_results;
extern std::vector<std::string> gene_names;

void internal_setup(std::string antimony_file);
Rcpp::NumericVector internal_evolution(Rcpp::NumericVector parent_state, Rcpp::Environment node);

#endif
