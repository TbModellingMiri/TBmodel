rm(list = ls())
library(ggplot2)
# shared functions ----
func_get_newist_prefix <- function(dir_name) {
    df <- file.info(list.files(output_dir, full.names = TRUE))
    newist_file <- rownames(df)[which.max(df$mtime)]
    newist_file <- tail(strsplit(newist_file, "/")[[1]], n = 1)
    return(head(strsplit(newist_file, "_")[[1]], n = 1))
}


func_drop_last_col <- function(df) {
    return(df[1:(length(df) - 1)])
}

##
# get newist run's prefix ----

output_dir <- file.path("ibm", "outputs")
newist_file_pf <- func_get_newist_prefix(output_dir)
newist_file_pf

## plot prevalence ----

df_prev <- read.csv(
    file = file.path(
        output_dir,
        paste0(newist_file_pf, "_prev.csv")
    ), header = FALSE
)
df_prev <- func_drop_last_col(df_prev)

df_prev_regions <- read.csv(
    file = file.path(
        output_dir,
        paste0(newist_file_pf, "_prev_regions.csv")
    )
)
df_prev_regions <- func_drop_last_col(df_prev_regions)
df_prev_regions[["days"]] <- 1:length(df_prev_regions[[1]])
df_prev_regions[["all"]] <- df_prev[[1]]

tail(df_prev_regions)
df_prev_regions_long <- tidyr::pivot_longer(
    df_prev_regions,
    names_to = "Villages",
    values_to = "Prevalence",
    -days
)
head(df_prev_regions_long)

plt_prev <- ggplot() +
    geom_line(
        data = df_prev_regions_long,
        mapping = aes(x = days, y = Prevalence, color = Villages)
    ) +
    geom_line(
        data = df_prev_regions,
        mapping = aes(x = days, y = all)
    )
# plt_prev

## plot SEIR counts ----

df_sim_table <- read.csv(
    file = file.path(
        output_dir,
        paste0(newist_file_pf, "_sim_table.csv")
    )
)
df_sim_table <- func_drop_last_col(df_sim_table)
head(df_sim_table)
df_sim_table_long <- tidyr::pivot_longer(
    df_sim_table,
    names_to = "Compartments",
    values_to = "Count",
    -day
)
head(df_sim_table_long)
plt_comp <- ggplot() +
    geom_line(
        data = df_sim_table_long,
        mapping = aes(x = day, y = Count, color = Compartments)
    )
# plt_comp

# plot all ----

library("patchwork")
plt_comp / plt_prev