library(tidyverse)

index_data <- tibble(id = 1:2500,
       date_of_diagnosis = sample(seq(as.Date("2015-01-01"), as.Date("2020-12-31"), by = "day"),
                                  size = 2500,
                                  replace = TRUE)
       )

contact_data <- tibble(id = sample(index_data$id, size = 10000, replace = TRUE),
       active_tb_first_screening = sample(c("TB", "LTBI","No TB"),
                                          size = 10000,
                                          prob = c(0.005, 0.005, 0.99),
                                          replace = TRUE),
       active_tb_second_screening = sample(c("TB", "LTBI","No TB"),
                                          size = 10000,
                                          prob = c(0.005, 0.005, 0.99),
                                          replace = TRUE),
       active_tb_third_screening = sample(c("TB", "LTBI","No TB"),
                                          size = 10000,
                                          prob = c(0.005, 0.005, 0.99),
                                          replace = TRUE),
       active_tb_fourth_screening = sample(c("TB", "LTBI","No TB"),
                                          size = 10000,
                                          prob = c(0.005, 0.005, 0.99),
                                          replace = TRUE),
       other_risk_factor = round(runif(10000, 1, 10))
       )


household = tibble(household = round(runif(342800,1, 60000)))

write_csv(household, "household.csv")
write_csv(index_data, "index_data.csv")
write_csv(contact_data, "contact_data.csv")