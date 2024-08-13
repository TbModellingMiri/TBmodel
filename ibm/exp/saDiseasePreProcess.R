# /////////////////////////////////////////////////
# Sensitivity analysis disease, pre process
# /////////////////////////////////////////////////

# This R script is used run the sensitivity analysis of the disease parameters.
# The script is currently configured to run 50 simulations per parameter
# combination, which may be adjusted by changing the value assigned to
# num_runs_per_setting

library(tidyverse)
library(doParallel)
library(foreach)
library(data.table)
library(rjson)

test_dir <- file.path("ibm/exp")
base_config <- rjson::fromJSON(file = file.path(test_dir, "base_config.json"))
exe <- file.path("ibm/exp", "ibm.exe")
config_file_name <- "config.json"

## Simulation setup
num_runs_per_setting <- 50
base_config$simulation$length_days <- 12652

## Parameter space
seq_hha <- c(0.02094982, 0.2, 0.4)
seq_beta <- c(0.1526, 0.1436, 0.1310)
seq_int <- c(1500, 1400, 1300)
seq_eta <- c(0.000005479452055*0.9,0.000005479452055*0.95, 0.000005479452055*1.05, 0.000005479452055*1.1)
seq_nu <- c(0.21*0.9,0.21*0.95, 0.21*1.05, 0.21*1.1)
seq_phi <- c(0.05*0.9,0.05*0.95, 0.05*1.05, 0.05*1.1)
seq_mu <- c(0.001065753*0.9,0.001065753*0.95, 0.001065753*1.05, 0.001065753*1.1)
seq_alpha <- c(0.00006191781*0.9,0.00006191781*0.95, 0.00006191781*1.05, 0.00006191781*1.1)
seq_tau <- c(573.05*0.9,573.05*0.95, 573.05*1.05, 573.05*1.1)

##### BEGIN #####
test_dir_outputs <- file.path(test_dir, "SA_phi")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

cc <- 0
for(i in seq_along(seq_hha)) {
    for(k in seq_along(seq_phi)) {
      
      config_new <- base_config
      cc <- cc + 1
      # Model parameters
      config_new$parameters$hh_attack_rate <- seq_hha[i]
      config_new$parameters$beta <- seq_beta[i]
      config_new$initial_condition$infectious$size <- seq_int[i]
      
      # SA
      config_new$parameters$phi <- seq_phi[k]
      
      
      config_new$data$path <- file.path(test_dir, config_new$data$path)
      
      config_path <- file.path(
        test_dir_outputs,
        paste0(
          "hha_", seq_hha[i], "_phi_", seq_phi[k]
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

##### END #####

##### BEGIN #####
test_dir_outputs <- file.path(test_dir, "SA_nu")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

cc <- 0
for(i in seq_along(seq_hha)) {
  for(k in seq_along(seq_nu)) {
    
    config_new <- base_config
    cc <- cc + 1
    # Model parameters
    config_new$parameters$hh_attack_rate <- seq_hha[i]
    config_new$parameters$beta <- seq_beta[i]
    config_new$initial_condition$infectious$size <- seq_int[i]
    
    # SA
    config_new$parameters$nu <- seq_nu[k]
    
    
    config_new$data$path <- file.path(test_dir, config_new$data$path)
    
    config_path <- file.path(
      test_dir_outputs,
      paste0(
        "hha_", seq_hha[i], "_nu_", seq_nu[k]
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

##### END #####

##### BEGIN #####
test_dir_outputs <- file.path(test_dir, "SA_eta")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

cc <- 0
for(i in seq_along(seq_hha)) {
  for(k in seq_along(seq_eta)) {
    
    config_new <- base_config
    cc <- cc + 1
    # Model parameters
    config_new$parameters$hh_attack_rate <- seq_hha[i]
    config_new$parameters$beta <- seq_beta[i]
    config_new$initial_condition$infectious$size <- seq_int[i]
    
    # SA
    config_new$parameters$eta <- seq_eta[k]
    
    
    config_new$data$path <- file.path(test_dir, config_new$data$path)
    
    config_path <- file.path(
      test_dir_outputs,
      paste0(
        "hha_", seq_hha[i], "_eta_", seq_eta[k]
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

##### END #####

##### BEGIN #####
test_dir_outputs <- file.path(test_dir, "SA_alpha")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

cc <- 0
for(i in seq_along(seq_hha)) {
  for(k in seq_along(seq_alpha)) {
    
    config_new <- base_config
    cc <- cc + 1
    # Model parameters
    config_new$parameters$hh_attack_rate <- seq_hha[i]
    config_new$parameters$beta <- seq_beta[i]
    config_new$initial_condition$infectious$size <- seq_int[i]
    
    # SA
    config_new$parameters$alpha <- seq_alpha[k]
    
    
    config_new$data$path <- file.path(test_dir, config_new$data$path)
    
    config_path <- file.path(
      test_dir_outputs,
      paste0(
        "hha_", seq_hha[i], "_alpha_", seq_alpha[k]
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

##### END #####

##### BEGIN #####
test_dir_outputs <- file.path(test_dir, "SA_mu")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

cc <- 0
for(i in seq_along(seq_hha)) {
  for(k in seq_along(seq_mu)) {
    
    config_new <- base_config
    cc <- cc + 1
    # Model parameters
    config_new$parameters$hh_attack_rate <- seq_hha[i]
    config_new$parameters$beta <- seq_beta[i]
    config_new$initial_condition$infectious$size <- seq_int[i]
    
    # SA
    config_new$parameters$mu <- seq_mu[k]
    
    
    config_new$data$path <- file.path(test_dir, config_new$data$path)
    
    config_path <- file.path(
      test_dir_outputs,
      paste0(
        "hha_", seq_hha[i], "_mu_", seq_mu[k]
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

##### END #####

##### BEGIN #####
test_dir_outputs <- file.path(test_dir, "SA_tau")

stopifnot(!dir.exists(test_dir_outputs))
dir.create(test_dir_outputs)

cc <- 0
for(i in seq_along(seq_hha)) {
  for(k in seq_along(seq_tau)) {
    
    config_new <- base_config
    cc <- cc + 1
    # Model parameters
    config_new$parameters$hh_attack_rate <- seq_hha[i]
    config_new$parameters$beta <- seq_beta[i]
    config_new$initial_condition$infectious$size <- seq_int[i]
    
    # SA
    config_new$parameters$tau <- seq_tau[k]
    
    
    config_new$data$path <- file.path(test_dir, config_new$data$path)
    
    config_path <- file.path(
      test_dir_outputs,
      paste0(
        "hha_", seq_hha[i], "_tau_", seq_tau[k]
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

##### END #####