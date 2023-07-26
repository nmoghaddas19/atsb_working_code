library(lubridate)
library(dplyr)
zambia_catch <- read.csv("~/Documents/GitHub/zambia_feeding_data")
zambia_catch$days <- (month(zambia_catch$collection_date)-1)*30 + day(zambia_catch$collection_date)

zambia_catch |>
  group_by(days) |>
  summarise(prop_fun = sum(mosquito_species == "funestus")/n(),
            prop_gamb = sum(mosquito_species == "gambiae")/n(),
            total_catch = n()) -> zambia_catch_clustered

bootstrap <- function(data) {
  statistic <- c()
  data_filtered <- zambia_catch_clustered
  size = length(data_filtered$prop_fun)
  for (i in 1:5000) {
    boot_sample <- sample(data_filtered$prop_fun, 
                          prob = data_filtered$total_catch, size = size,
                          replace = TRUE)
    statistic[i] <- mean(boot_sample)
  }
  return(
    list(
      statistic = statistic,
      quantiles = quantile(statistic, probs = c(0.025,0.975)),
      size = size,
      mean = mean(statistic)
    )
  )
}
vector_proportions <- bootstrap(zambia_catch_clustered)









# zambia_catch |>
#   group_by(cluster) |>
#   summarise(prop_fun = sum(mosquito_species == "funestus")/n(),
#             prop_gamb = sum(mosquito_species == "gambiae")/n()) -> zambia_catch_clustered
# 
# bootstrap <- function(data, clust) {
#   statistic <- c()
#   data |>
#     filter(cluster == clust) |>
#     group_by(days) |>
#     summarise(prop_fun = sum(mosquito_species == "funestus")/n(),
#               prop_gamb = sum(mosquito_species == "gambiae")/n(),
#               total_catch = n()) -> data_filtered
#   size = length(data_filtered$prop_fun)
#   for (i in 1:5000) {
#     boot_sample <- sample(data_filtered$prop_fun, 
#                           prob = data_filtered$total_catch, size = size,
#                           replace = TRUE)
#     statistic[i] <- mean(boot_sample)
#   }
#   return(
#     list(
#       statistic = statistic,
#       quantiles = quantile(statistic, probs = c(0.025,0.975)),
#       size = size,
#       cluster = clust,
#       mean = mean(statistic)
#     )
#   )
# }
# 
# clusters <- unique(zambia_catch$cluster)
# vector_proportions <- list()
# for (i in 1:length(clusters)) {
#   vector_proportions[[clusters[i]]] <- bootstrap(zambia_catch, clusters[i])
# }