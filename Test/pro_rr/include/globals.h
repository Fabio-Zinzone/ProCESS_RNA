#ifndef GLOBALS_H
#define GLOBALS_H

#include <cstdint>
#include <string>
#include <vector>
#include <Rcpp.h>
#include <rr/rrRoadRunner.h>
#include <rr/rrExecutableModel.h>
#include "antimony_api.h"

using namespace rr;
using namespace std;
using namespace Rcpp;

extern RoadRunner* rrg;

void internal_setup(std::string antimony_file);
Rcpp::NumericVector internal_evolution(Rcpp::NumericVector parent_state, Rcpp::Environment node);

#endif
