# /////////////////////////////////////////////////
# Fit, pre process
# /////////////////////////////////////////////////

# This R script is used run fitting proceedure. The script is configured to
# run 40 simulations per parameter combination, which may be adjusted by
# changing the value assigned to num_runs_per_setting.
# The granularity and bounds of the grid search can by adjusting the beta_seqs

library(tidyverse)
library(doParallel)
library(foreach)
library(data.table)
library(rjson)

test_dir <- file.path("ibm/exp")
base_config <- rjson::fromJSON(file = file.path(test_dir, "base_config.json"))
exe <- file.path("ibm/exp/ibm.exe")
config_file_name <- "config.json"

## Simulation setup
num_runs_per_setting <- 40
base_config$simulation$length_days <- 9000

## Parameter space
beta_seqs <- list()
beta_seqs[[1]] <- seq(from = 0.125, to = 0.165, by = 0.005)
beta_seqs[[2]] <- seq(from = 0.125, to = 0.165, by = 0.005)
beta_seqs[[3]] <- seq(from = 0.125, to = 0.165, by = 0.005)
seq_hha <- c(0.02094982, 0.2, 0.4)
seq_int <- c(1500, 1400, 1300)

## Create space
test_dir_outputs <- file.path(test_dir, "fit")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

cc <- 0
for(ihha in seq_along(seq_hha)) {
  seq_beta <- beta_seqs[[ihha]]
  for(ibeta in seq_along(seq_beta)) {
    config_new <- base_config
    cc <- cc + 1                
    config_new$parameters$hh_attack_rate <- seq_hha[ihha]
    config_new$initial_condition$infectious$size <- seq_int[ihha]
    
    config_new$parameters$beta <- seq_beta[ibeta]
    
    config_new$data$path <- file.path(test_dir, config_new$data$path)
    
    config_path <- file.path(
      test_dir_outputs,
      paste0(
        "hha_", seq_hha[ihha], "_beta_", seq_beta[ibeta]
      ))
    
    dir.create(config_path)
    
    write(
      rjson::toJSON(config_new),
      file = file.path(
        config_path,
        config_file_name
      )
    )
  }
}



############## RUNNER ##############
num_cores <- detectCores() - 1 # !! doParallel::detectCores() detects the number of cores. 
                               # on a laptop, this detects the number of cores on a single CPU
                               # I don't know how this detects core across nodes on a HPC
registerDoParallel(num_cores)

time_start_parallel <- Sys.time()

all_configs_dirs <- list.dirs(test_dir_outputs, recursive = FALSE)

foreach(iic = seq_along(all_configs_dirs), .combine='c') %:%
  foreach(iir = 1:num_runs_per_setting, .combine='c') %dopar% {
    
    cmd <- paste(
      exe,
      "-c", file.path(all_configs_dirs[iic], config_file_name),
      "-o", file.path(all_configs_dirs[iic], "outputs"),
      "-p", paste0("c_", iic, "r_", iir)
    )
    system(cmd)
  }

# Stop parallel backend
time_stop_parallel <- Sys.time()
print(time_stop_parallel-time_start_parallel)
stopImplicitCluster() 