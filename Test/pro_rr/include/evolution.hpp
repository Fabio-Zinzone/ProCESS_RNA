#include "globals.h"

NumericVector rr_evolution(NumericVector parent_state, Environment node) {
  if (rrg == nullptr) stop("Roadrunner not initialized");
  
  double duration = (double)node["death_time"] - (double)node["birth_time"];
  if (duration <= 0.0) return parent_state;
  
  rrg->reset();
  CharacterVector names = parent_state.names();
  
  for (int i = 0; i < parent_state.size(); ++i) {
    rrg->setValue(as<string>(names[i]), parent_state[i]);
  }
  
  const ls::DoubleMatrix* result = rrg->simulate(0.0, duration, 2);
  
  int last_row = result->RSize() - 1;
  NumericVector final_state(result->CSize());
  
  for (int j = 0; j < result->CSize(); ++j) {
    final_state[j] = (*result)(last_row, j);
  }
  
  final_state.names() = names;
  return final_state;
}
