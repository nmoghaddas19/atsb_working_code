library(netz)
library(dplyr)
library(foresite)
library(site)
library(malariasimulation)

cluster_summaries <- read.csv("~/Documents/GitHub/atsb_working_code/Zambia ATSB projections/baseline_cluster_summaries_9.7.21 (for restricted randomization n=70).csv")
pyr_nets <- read.csv("~/Documents/GitHub/pyrethroid_only_nets.csv")
pbo_nets <- read.csv("~/Documents/GitHub/pyrethroid_pbo_nets.csv")

data <- cluster_summaries[,c(1,2,3,6,7)]
data$baseline_date <- rep("09.07.2021", 70)
dist_baseline_interval <- as.numeric(abs(difftime(dmy("15 Jan 2021"), 
                                                  dmy(data$baseline_date)[1], 
                                                  units = "days")))

distribution_tt <- 365 + 15
target_tt <- 365 + 15 + dist_baseline_interval
jan_coverages <- numeric(70)
for (i in 1:nrow(data)) {
  target <- c(data$net_use[i]/100)
  fit <- fit_usage(target_usage = target, target_usage_timesteps = target_tt,
                   distribution_timesteps = distribution_tt)
  jan_coverages[i] <- fit$par
}
data$jan2021_coverage <- jan_coverages
data$mar2022_coverage <- jan_coverages # no data (yet) to base this off so assume same as previous distribution
data$oct2022_coverage <- jan_coverages # no data (yet) to base this off so assume same as previous distribution

itn_distributions <- data.frame(
  date = c("January 2021", "March 2022", "October 2022"),
  type = c("Pyrethroid", "Pyrethroid", "PBO")
)
irs_campaigns <- data.frame(
  date = c("December 2021", "December 2022"),
  type = c("Fludora fusion", "Fludora fusion")
)
clothianidin <- data.frame(
  ls_theta = c(0.792, -0.127, 2.035),
  ls_gamma = c(-0.007, -0.010, -0.006), 
  ks_theta = c(-1.328, -2.481, -0.556),
  ks_gamma = c(0.009, 0.007, 0.011),
  ms_theta = c(-0.458, -1.445, -0.092),
  ms_gamma = c(-0.007, -0.010, -0.006) 
)
rownames(clothianidin) <- c("mean", "lower", "upper")

seasonality <- single_site(ZMB, 18)$seasonality

for (i in 1:70) {
  cluster_params <- get_parameters(list(
    human_population = 10000,
    model_seasonality = TRUE,
    g0 = seasonality$g0,
    g = c(seasonality$g1, seasonality$g2, seasonality$g3),
    h = c(seasonality$h1, seasonality$h2, seasonality$h3),
    individual_mosquitoes = FALSE,
    atsb = FALSE
  ))
  index <- round(runif(1, min=1, max=5000))
  cluster_params <- set_species(parameters = cluster_params,
                                species = list(fun_params, gamb_params),
                                proportions = c(vector_proportions$statistic[index], 
                                                1-vector_proportions$statistic[index]))
  fun_res <- 100-round(funestus_resistance$statistic[index])
  gam_res <- 100-round(gambiae_resistance$statistic[index])
  dn0_pyr_fun <- pyr_nets$dn0_med[fun_res+1]
  dn0_pbo_fun <- pbo_nets$dn0_med[fun_res+1]
  rn0_pyr_fun <- pyr_nets$rn0_med[fun_res+1]
  rn0_pbo_fun <- pbo_nets$rn0_med[fun_res+1]
  gamman_pyr_fun <- pyr_nets$gamman_med[fun_res+1]
  gamman_pbo_fun <- pbo_nets$gamman_med[fun_res+1]
  dn0_pyr_gam <- pyr_nets$dn0_med[gam_res+1]
  dn0_pbo_gam <- pbo_nets$dn0_med[gam_res+1]
  rn0_pyr_gam <- pyr_nets$rn0_med[gam_res+1]
  rn0_pbo_gam <- pbo_nets$rn0_med[gam_res+1]
  gamman_pyr_gam <- pyr_nets$gamman_med[gam_res+1]
  gamman_pbo_gam <- pbo_nets$gamman_med[gam_res+1]
}
