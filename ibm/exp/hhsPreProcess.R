# /////////////////////////////////////////////////
# Household screening, pre process
# /////////////////////////////////////////////////

# This R script is used run the simulation of the household screening
# intervention. The script is currently configured to run 50 simulations per
# parameter combination, which may be adjusted by changing the value assigned
# to num_runs_per_setting

library(tidyverse)
library(doParallel)
library(foreach)
library(data.table)
library(rjson)

test_dir <- file.path("ibm/exp")
base_config <- rjson::fromJSON(file = file.path(test_dir, "base_config.json"))
exe <- file.path("ibm/exp", "ibm.exe")
config_file_name <- "config.json"
test_dir_outputs <- file.path(test_dir, "intervention_hhs")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

## Simulation seutp
num_runs_per_setting <- 50
base_config$simulation$length_days <- 12652

## parameter range
seq_beta <- c(0.1526, 0.1436, 0.1310)
seq_hha <- c(0.02094982, 0.2, 0.4)
seq_int <- c(1500, 1400, 1300)
seq_hhs <- c(0, 1)



## Create space
cc <- 0
for (i in seq_along(seq_hha)) {
  for(j in seq_along(seq_hhs)) {
    
    config_new <- base_config
    cc <- cc + 1
    # Model parameters
    config_new$parameters$hh_attack_rate <- seq_hha[i]
    config_new$parameters$beta <- seq_beta[i]
    config_new$initial_condition$infectious$size <- seq_int[i]
    
    # Internvetions
    config_new$intervention$HHS <- seq_hhs[j]

    config_new$data$path <- file.path(test_dir, config_new$data$path)
    
    config_path <- file.path(
      test_dir_outputs,
      paste0(
        "hha_", seq_hha[i], "_hhs_", seq_hhs[j]
      )
    )
    
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
num_cores <- detectCores() - 1
registerDoParallel(num_cores)

time_start_parallel <- Sys.time()

all_configs_dirs <- list.dirs(test_dir_outputs, recursive = FALSE)

foreach(iic = seq_along(all_configs_dirs), .combine = "c") %:%
  foreach(iir = 1:num_runs_per_setting, .combine = "c") %dopar% {
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
print(time_stop_parallel - time_start_parallel)
stopImplicitCluster()

