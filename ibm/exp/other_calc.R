# /////////////////////////////////////////////////
# Other calculations
# /////////////////////////////////////////////////

# This R script is used make the supplementary calculations for the model.

library(tidyverse)

index_data <- read_csv("ibm/exp/data/index_data.csv")
contact_data <- read_csv("ibm/exp/data/contact_data.csv")

# calculate SAR
index_data %>% 
  left_join(contact_data) %>% 
  mutate(any_active = apply(select(., contains("active")), 1,
                            function(.x) {any(.x %in% c("TB","LTBI"), na.rm = TRUE)}), 
         household = ifelse(other_risk_factor %in% c(1,2,3), #  is code for household contact
                            TRUE, FALSE)) %>% 
  filter(household == TRUE) %>% 
  reframe(sar = mean(any_active) , .by = id) %>% 
  reframe(sar = mean(sar))

# calculate household distribution
index_data %>% 
  left_join(contact_data) %>% 
  count(id) %>% 
  mutate(n = n+1) %>% 
  pull()
