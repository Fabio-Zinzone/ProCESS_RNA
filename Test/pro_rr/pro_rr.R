suppressPackageStartupMessages({
library(Rcpp)
library(dplyr)
library(ProCESS)
})

message("SIMULATION STARTED")

start_time <- Sys.time()

sourceCpp("master.cpp")
set.seed(0)
set_cpp_seed(0)

sim <- TissueSimulation()

sim$add_mutant(name = "A", growth_rates = 0.5, death_rates = 0.01)

sim$place_cell("A", 500, 500)

sim$run_up_to_size("A", 1000)

plot_tissue(sim)

sim$sample_cells("S", num_of_cells = 10)

forest <- sim$get_sample_forest()

plot_forest(forest)

roadrunner_setup("grn_anti_10.txt")

initial_grn_state <- c(G0 = 5.0, G1 = 1.0, G2 = 0.5, G3 = 0.5, G4 = 2.0, 
                       G5 = 5.0, G6 = 1.0, G7 = 0.5, G8 = 4.0, G9 = 0.2)

tour <- get_label_tour(forest, 
                      roadrunner_evolution, 
                      init_value = initial_grn_state, 
                      only_leaves = FALSE)

while (!tour$done) {
  tour$step()
}

apply_poisson()
results <- save_grn_matrix()
#storage.mode(results)
save_csv_matrix("count_matrix.csv")
message("Count Matrix saved as CSV")

end_time <- Sys.time()
delta_time <- end_time - start_time
print(delta_time)