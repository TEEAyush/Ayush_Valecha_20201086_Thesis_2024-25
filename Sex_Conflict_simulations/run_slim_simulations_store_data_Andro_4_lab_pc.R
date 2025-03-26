# # Load required packages
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages("ggplot2")




library(doParallel)
library(foreach)
library(dplyr)
library(stringr)
library(ggplot2)


source("/mnt/data2/avalecha/slim/SLiM/My_scripts/Andro_4/functions_Andro_4_lab_pc.R")


#################################################
# s and r values
s_values <- c(0, 0.1, 0.2, 0.3, 0.5)
#s_values <- c(0.1, 0.5)
#s_values <- c(0)
r_values <- c(0.001, 0.01, 0.1, 0.5)
#r_values <- c(0.001, 0.5)
#r_values <- c(0.5)
#num_reps <- 2
num_reps <- 10
###################################################




# Create a grid of all combinations
param_grid <- expand.grid(
  s = s_values,
  r = r_values,
  rep = 1:num_reps,
  stringsAsFactors = FALSE
)



# Define the directory to store simulation results
#results_dir <- "mnt/data2/avalecha/slim/SLiM/Test_Simulation_Results/"
results_dir <- "/mnt/data2/avalecha/slim/SLiM/Simulation_Results_Andro_4/beta_2_20_reps_2/neutral_results_full/"

# Create the directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

###################################################
# Other parameters
L <- 100
#t <- 500
mut_init_freq <- 500
K <- 500  # or 1000, up to you
##########################################################sim


# Detect the number of available cores
num_cores <- parallel::detectCores()
cat("Number of available cores:", num_cores, "\n")

# Decide how many cores to use (leave one free)
cores_to_use <- 8
cat("Using", cores_to_use, "cores for parallel processing.\n")

# Register the parallel backend
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)

# Start parallel simulations
results_list <- foreach(i = 1:nrow(param_grid), .packages = c("dplyr", "stringr"), .export = c("run_slim_simulation2", "parse_slim2_output")) %dopar% {
  
  # Extract parameter values for this iteration
  current_s <- param_grid$s[i]
  current_r <- param_grid$r[i]
  current_rep <- param_grid$rep[i]
  
  # Define a unique filename based on parameters
  filename <- paste0(
    "neutral_sim_s_", current_s, 
    "_r_", current_r, 
    "_rep_", current_rep,
    ".rds"
  )
  
  
  # Full path for the output file
  file_path <- file.path(results_dir, filename)
  
  # Run the simulation
  sim_result <- run_slim_simulation2(
    s = current_s,
    r = current_r,
    K = K,
    L = L,
    mut_init_freq = mut_init_freq
  )
  
  # Save the simulation result to disk
  saveRDS(sim_result, file = file_path)
  
  # Optionally, return NULL to reduce memory usage
  return(NULL)
}


# Stop the cluster after simulations are done
stopCluster(cl)

cat("All simulations completed and results saved to disk.\n")