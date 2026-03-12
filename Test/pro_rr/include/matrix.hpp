#include "globals.h"

std::map<int, std::vector<double>> grn_results;
std::vector<std::string> gene_names;

int poisson(double lambda) {
  if (lambda <= 0.0) return 0;
  std::poisson_distribution<int> d(lambda);
  return d(gen);
}

void int_poisson() {
  for (auto& pair : grn_results) {
    for (double& val : pair.second) {
      val = poisson(val);
    }
  }
}

Matrix_type grn_matrix() {
  
  if (grn_results.empty()) { Rcpp::stop("No data available"); }
  
  int n_cells = grn_results.size();
  int n_genes = gene_names.size();
  
  Matrix_type mat(n_genes, n_cells);
  Rcpp::CharacterVector col_names(n_cells);
  Rcpp::CharacterVector row_names(n_genes);
  
  for (int r = 0; r < n_genes; ++r) {
      row_names[r] = gene_names[r];
  }
  
  int col_idx = 0;
  for (auto const& pair : grn_results) {
    int cell_id = pair.first;
    const std::vector<double>& genes = pair.second;
     
    col_names[col_idx] = std::to_string(cell_id);
      
    for (int row_idx = 0; row_idx < n_genes; ++row_idx) {
      mat(row_idx, col_idx) = genes[row_idx];
    }
    col_idx++;
  }
  
  mat.attr("dimnames") = Rcpp::List::create(row_names, col_names);
  
  return mat;
}

void save_matrix(std::string filename) {
  
  if (grn_results.empty()) { Rcpp::stop("No data available"); }
  
  std::ofstream file(filename);
  
  if (!file.is_open()) { Rcpp::stop("Error: impossible csv creation"); }
  
  int n_genes = gene_names.size();
  
  file << " ";
  for (auto const& pair : grn_results) {
    file << "," << pair.first;
  }
  file << "\n";
  
  for (int r = 0; r < n_genes; ++r) {
    file << gene_names[r];
        
    for (auto const& pair : grn_results) {
        file << "," << pair.second[r];
    }
    file << "\n";
  }
  
  file.close();
}
