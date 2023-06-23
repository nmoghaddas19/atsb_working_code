# This script will draw heavily on my ASB workflow to generate estimates 
# of the uncertainty in prevalence caused by not measuring feeding rates
# in trial clusters. I will do this by assuming the variability in dyed 
# fraction in control clusters is comparable to what would be seen in the
# trial clusters. 

library(malariasimulation)
library(dplyr)
library(RColorBrewer)
library(foresite)
library(site)
library(lubridate)

# read in the dyed fraction data
zambia <- read.csv("~/Documents/GitHub/atsb_working_code/zambia_asb_data")
mali <- read.csv("~/Documents/GitHub/atsb_working_code/mali_asb_data")

# extract static cluster and species specific dyed fractions
# mali has no data on species specific dye positivity therefore we will pool
# the data
mali |>
  group_by(Village) |>
  summarise(dyed_fraction=weighted.mean(dyed_fraction,total_sampled),
            total_catch=sum(total_catch)) -> mali_grouped
villages <- mali_grouped$Village
cluster_quantiles <- matrix(0,nrow = length(villages), ncol=2)
for (i in 1:length(villages)) {
  this_village <- villages[i]
  mali_filtered <- filter(mali, Village==this_village)
  replicates <- matrix(0,nrow=5000,ncol=nrow(mali_filtered))
  for (j in 1:5000) {
    replicates[j,] <- sample(x = mali_filtered$total_asb_positive/mali_filtered$total_sampled, 
                             replace = TRUE, 
                             size = nrow(mali_filtered), 
                             prob = mali_filtered$total_sampled)
  }
  cluster_means <- apply(replicates, MARGIN = 1, FUN = mean)
  cluster_quantiles[i,] <- quantile(cluster_means, probs = c(0.025,0.975))
}
colnames(cluster_quantiles) <- c("feed_rate_lower", "feed_rate_upper")
mali_grouped <- cbind(mali_grouped, cluster_quantiles)
par(las=2, mfrow=c(1,1), mar=c(8,4,4,1)+0.1)
plot(mali_grouped$dyed_fraction,
     col=brewer.pal(7,"Set2")[factor(mali_grouped$Village)],
     pch=19,
     cex=2*mali_grouped$total_catch/mean(mali_grouped$total_catch),
     frame.plot = F,
     ylim = c(0,1),
     ylab = "Dyed fraction",
     xaxt = "n",
     xlab = "")
axis(1, at=1:7, labels=mali_grouped$Village)
grid()
arrows(1:7,
       mali_grouped[,4],
       1:7,
       mali_grouped[,5],
       code=0,
       lwd=2,
       col = brewer.pal(7,"Set2")[factor(mali_grouped$Village)])

zambia$days <- (month(zambia$collection_date)-1)*30 + day(zambia$collection_date)

zambia |>
  group_by(cluster, mosquito_species, days) |>
  summarise(total_catch = n(),
            dyed_fraction = mean(positive)) -> zambia_grouped

villages <- unique(zambia_grouped$cluster)
cluster_quantiles <- matrix(0,nrow = length(villages), ncol=2)
for (i in 1:length(villages)) {
  this_village <- villages[i]
  zambia_filtered <- filter(zambia_grouped, cluster==this_village, mosquito_species == "funestus")
  replicates <- matrix(0,nrow=5000,ncol=nrow(zambia_filtered))
  for (j in 1:5000) {
    replicates[j,] <- sample(x = zambia_filtered$dyed_fraction, 
                             replace = TRUE, 
                             size = nrow(zambia_filtered), 
                             prob = zambia_filtered$total_catch)
  }
  cluster_means <- apply(replicates, MARGIN = 1, FUN = mean)
  cluster_quantiles[i,] <- quantile(cluster_means, probs = c(0.025,0.975))
}
colnames(cluster_quantiles) <- c("feed_rate_lower_f", "feed_rate_upper_f")

cluster_quantiles_gambiae <- matrix(0,nrow = length(villages), ncol=2)
for (i in 1:length(villages)) {
  this_village <- villages[i]
  zambia_filtered <- filter(zambia_grouped, cluster==this_village, mosquito_species == "gambiae")
  replicates <- matrix(0,nrow=5000,ncol=nrow(zambia_filtered))
  for (j in 1:5000) {
    replicates[j,] <- sample(x = zambia_filtered$dyed_fraction, 
                             replace = TRUE, 
                             size = nrow(zambia_filtered), 
                             prob = zambia_filtered$total_catch)
  }
  cluster_means <- apply(replicates, MARGIN = 1, FUN = mean)
  cluster_quantiles_gambiae[i,] <- quantile(cluster_means, probs = c(0.025,0.975))
}
colnames(cluster_quantiles_gambiae) <- c("feed_rate_lower_g", "feed_rate_upper_g")

zambia |>
  group_by(cluster, mosquito_species) |>
  summarise(total_catch = n(),
            dyed_fraction = mean(positive),
            .groups = ) -> zambia_grouped
plot(
  zambia_grouped$dyed_fraction,
  col=brewer.pal(7,"Set2")[factor(zambia_grouped$mosquito_species)],
  pch=19,
  cex=1.2,
  frame.plot = F,
  ylim = c(0,0.5),
  ylab = "Dyed fraction",
  xaxt = "n",
  xlab = ""
)
grid()
# # Function to calculate the weighted mean by cluster
# weighted_mean_by_cluster <- function(data) {
#   cluster_weights <- tapply(data$total_catch, data$cluster, sum)
#   cluster_means <- tapply(data$dyed_fraction * data$total_catch, data$cluster, sum) / cluster_weights
#   return(cluster_means)
# }
# 
# # Function to calculate the bootstrapped mean by cluster
# bootstrapped_mean_by_cluster <- function(data) {
#   bootstrapped_means <- tapply(data$dyed_fraction, data$cluster, function(x) {
#     indices <- sample(length(x), replace = TRUE)
#     return(mean(x[indices]))
#   })
#   return(bootstrapped_means)
# }
# 
# # Number of bootstrap iterations
# n_iterations <- 1000
# 
# zambia_grouped |> filter(mosquito_species=="funestus") ->cluster_data
# # Calculate weighted mean by cluster
# cluster_weighted_means <- weighted_mean_by_cluster(cluster_data)
# 
# # Calculate bootstrapped means by cluster
# cluster_bootstrapped_means <- bootstrapped_mean_by_cluster(cluster_data)
# 
# # Function to calculate bootstrap confidence intervals
# bootstrap_ci <- function(data) {
#   lower_ci <- quantile(data, 0.025)
#   upper_ci <- quantile(data, 0.975)
#   return(c(lower_ci, upper_ci))
# }
# 
# # Calculate bootstrap confidence intervals by cluster
# cluster_ci <- tapply(cluster_bootstrapped_means, names(cluster_bootstrapped_means), bootstrap_ci)
