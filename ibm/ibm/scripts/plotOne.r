library(tidyverse)
library(doParallel)
library(foreach)
my_dir <- file.path("~/Miri/ibm/ibm")
setwd(my_dir)

exe <- file.path("ibm.exe")
config_file_name <- file.path("config.json")
output_folder <- file.path("plot_one")
name <- "plot"


stopifnot(!dir.exists(output_folder))
dir.create(output_folder)

num_cores <- detectCores() - 1
registerDoParallel(num_cores)

runs<-num_cores

{ 
start <- Sys.time()
foreach(i = 1:runs, .combine='c') %dopar% {
  cmd <- paste(
    exe,
    "-c", file.path(config_file_name),
    "-o", file.path(output_folder),
    "-p", i
  )
  system(cmd)
}
stop <- Sys.time()
elapsed <- stop-start
print(elapsed)



csv <- list.files(output_folder, pattern = "sim_table.csv")
csvs <- list()

for(i in seq_along(csv)) {
  df <- as_tibble(data.table::fread(paste0(output_folder,"/", csv[i])))
  df <- df[, -ncol(df)]
  df$index <- i
  csvs[[i]] <- df
}

data <- as_tibble(do.call(rbind, csvs))

run_in <- 6000
sim_length <- 9000

sim_start_date <- as.Date("2021-01-01")-run_in

inc_data <- data %>% 
  mutate(index = (row_number()-1) %/% (sim_length-1)) %>% 
  mutate(date = sim_start_date+day,
         year = substr(date, 1,4),
         month = paste0(substr(date, 6,7),"-","01"),
         year_month = as.Date(paste0(year, "-",month))) %>%  
  reframe(inc = sum(inc),.by = c(year_month, index))

mean_inc_data <- inc_data %>% 
  reframe(inc = mean(inc), .by = year_month)

inc_data %>%
  ggplot() +
  geom_line(data = inc_data,
            aes(year_month,
                inc, group = index), alpha = 0.1, color = "grey") +
  geom_line(data = mean_inc_data,
            aes(year_month,
                inc, col = "Mean")) + 
  labs(y = "Monthly incident",
       x = "Date",
       title = "Simulated effect of perfect household tracing",
       subtitle = expression(paste(beta, "=0.15; household attack rate = 0.5"))) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  scale_color_manual(values = c("red", "black"))

}
