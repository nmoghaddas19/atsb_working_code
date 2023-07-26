library(foresite)
library(site)
library(malariasimulation)
library(dplyr)
library(grr)

ig2_nets <- read.csv("~/Documents/GitHub/pyrethroid_pyrrole_nets.csv")

countries <- list(UGA, GHA, NGA)

names(countries) <- c("UGA", "GHA", "NGA")

uganda_districts <- data.frame(
  district = c("Bugiri", "Tororo", "Namutumba", "Butaleja", "Kibuku", "Budaka", 
               "Butebo", "Pallisa", "Serere", "Amolatar", "Kaberamaido", 
               "Kalaki", "Dokolo", "Alebtong", "Otuke", "Lira"),
  admin_1 = c("Bugiri", "Tororo", "Iganga", "Tororo", "Pallisa", "Pallisa",
              "Pallisa", "Pallisa", "Soroti", "Lira", "Kaberamaido", 
              "Kaberamaido", "Lira", "Lira", "Lira", "Lira"))

ghana_districts <- data.frame(
  district = c("Builsa North", "Builsa South", "Kasina Nankana West", "Daffiama Bussie",
               "Nadowli", "Jirapa", "Lawra", "Lambussie Karni.", "Nandom.", 
               "Sissala East.", "Sissala West", "Wa East", "Wa Municipal", 
               "Wa West", "Obuasi West", "Obuasi East", "West Mamprusi Districts", 
               "Gushegeu", "Karaga.", "Kumbungu.", "Mamprugu Moagduri.", "Bunkpurugu", 
               "Yunyoo", "Tatale", "East Mamprusi.", "Chereponi"),
  admin_1 = c("Upper East", "Upper East", "Upper East", "Upper West", "Upper West",
             "Upper West", "Upper West", "Upper West", "Upper West", "Upper West",
             "Upper West", "Upper West", "Upper West", "Upper West", "Ashanti",
             "Ashanti", "Northern", "Northern", "Northern", "Northern", "Northern",
             "Northern", "Northern", "Northern", "Northern", "Northern"))

foo <- matches(uganda_districts$admin_1, UGA$sites$name_1, all.y = FALSE)
foo <- foo[order(foo$x),]

data <- list()
for (i in 1:length(countries)) {
  country_site_files <- countries[[names(countries)[i]]]
  data[[names(countries)[i]]] <- run_country(country_site_files)
}
run_country(UGA)

run_country <- function(country_site_files) {
  country_name <- country_site_files$country
  country_results <- list()
  for (j in 1:length(foo$y)) {
    site <- single_site(country_site_files, foo$y[j])
    params <- site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir = site$eir$eir[1],
      overrides = list(human_population = 5000)
    )
    params$bednet_timesteps[24] <- params$bednet_timesteps[23] + 365
    params$bednet_coverages[24] <- 0.6
    params$bednet_dn0 <- rbind(params$bednet_dn0, params$bednet_dn0[23,])
    params$bednet_rn <- rbind(params$bednet_rn, params$bednet_rn[23,])
    params$bednet_rnm <- rbind(params$bednet_rnm, params$bednet_rnm[23,])
    params$bednet_gamman[24] <- params$bednet_gamman[23]
    
    params$spraying_timesteps[24] <- params$spraying_timesteps[23] + 365
    params$spraying_coverages[24] <- params$spraying_coverages[23]
    params$spraying_ls_theta <- rbind(params$spraying_ls_theta, params$spraying_ls_theta[23,])
    params$spraying_ls_gamma <- rbind(params$spraying_ls_gamma, params$spraying_ls_gamma[23,])
    params$spraying_ks_theta <- rbind(params$spraying_ks_theta, params$spraying_ks_theta[23,])
    params$spraying_ks_gamma <- rbind(params$spraying_ks_gamma, params$spraying_ks_gamma[23,])
    params$spraying_ms_theta <- rbind(params$spraying_ms_theta, params$spraying_ms_theta[23,])
    params$spraying_ms_gamma <- rbind(params$spraying_ms_gamma, params$spraying_ms_gamma[23,])
    
    set.seed(j)
    control <- run_simulation(timesteps = params$timesteps + 5*365, 
                              parameters = params)
    
    feeding_rate <- runif(1, 0.05, 0.35)
    params_atsb_noirs <- site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir = site$eir$eir[1],
      overrides = list(human_population = 5000,
                       atsb = TRUE,
                       mu_atsb = c(feeding_rate, feeding_rate, feeding_rate))
    )
    params_atsb_noirs <- set_atsb(
      parameters = params_atsb_noirs,
      timesteps = (params_atsb_noirs$spraying_timesteps[23] + 365):(params_atsb_noirs$spraying_timesteps[23] + 2*365), 
      coverages = rep(1,366)
      )
    
    params_atsb_noirs$bednet_timesteps[24] <- params_atsb_noirs$bednet_timesteps[23] + 365
    params_atsb_noirs$bednet_coverages[24] <- 0.6
    params_atsb_noirs$bednet_dn0 <- rbind(params_atsb_noirs$bednet_dn0, params_atsb_noirs$bednet_dn0[23,])
    params_atsb_noirs$bednet_rn <- rbind(params_atsb_noirs$bednet_rn, params_atsb_noirs$bednet_rn[23,])
    params_atsb_noirs$bednet_rnm <- rbind(params_atsb_noirs$bednet_rnm, params_atsb_noirs$bednet_rnm[23,])
    params_atsb_noirs$bednet_gamman[24] <- params_atsb_noirs$bednet_gamman[23]
    
    set.seed(j)
    atsb_noirs <- run_simulation(timesteps = params_atsb_noirs$timesteps + 5*365, 
                                 parameters = params_atsb_noirs)
    
    params_bells_whistles <- site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir = site$eir$eir[1],
      overrides = list(human_population = 5000,
                       atsb = TRUE,
                       mu_atsb = c(feeding_rate, feeding_rate, feeding_rate))
    )
    params_bells_whistles <- set_atsb(
      parameters = params_bells_whistles,
      timesteps = (params_bells_whistles$spraying_timesteps[23] + 365):(params_bells_whistles$spraying_timesteps[23] + 2*365), 
      coverages = rep(1,366)
    )
    res <- round(0)
    params_bells_whistles$bednet_timesteps[24] <- params_bells_whistles$bednet_timesteps[23] + 365
    params_bells_whistles$bednet_coverages[24] <- 0.6
    params_bells_whistles$bednet_dn0 <- rbind(params_bells_whistles$bednet_dn0, rep(ig2_nets$dn0_med[1],3))
    params_bells_whistles$bednet_rn <- rbind(params_bells_whistles$bednet_rn, rep(ig2_nets$rn0_med[1],3))
    params_bells_whistles$bednet_rnm <- rbind(params_bells_whistles$bednet_rnm, params_bells_whistles$bednet_rnm[23,])
    params_bells_whistles$bednet_gamman[24] <- ig2_nets$gamman_med[1]*365
    
    params_bells_whistles$spraying_timesteps[24] <- params_bells_whistles$spraying_timesteps[23] + 365
    params_bells_whistles$spraying_coverages[24] <- params_bells_whistles$spraying_coverages[23]
    params_bells_whistles$spraying_ls_theta <- rbind(params_bells_whistles$spraying_ls_theta, params_bells_whistles$spraying_ls_theta[23,])
    params_bells_whistles$spraying_ls_gamma <- rbind(params_bells_whistles$spraying_ls_gamma, params_bells_whistles$spraying_ls_gamma[23,])
    params_bells_whistles$spraying_ks_theta <- rbind(params_bells_whistles$spraying_ks_theta, params_bells_whistles$spraying_ks_theta[23,])
    params_bells_whistles$spraying_ks_gamma <- rbind(params_bells_whistles$spraying_ks_gamma, params_bells_whistles$spraying_ks_gamma[23,])
    params_bells_whistles$spraying_ms_theta <- rbind(params_bells_whistles$spraying_ms_theta, params_bells_whistles$spraying_ms_theta[23,])
    params_bells_whistles$spraying_ms_gamma <- rbind(params_bells_whistles$spraying_ms_gamma, params_bells_whistles$spraying_ms_gamma[23,])
    
    set.seed(j)
    atsb_bells_whistles <- run_simulation(timesteps = params_bells_whistles$timesteps + 5*365, 
                                          parameters = params_bells_whistles)
    site_name <- paste0(site$sites$name_1, "_", site$sites$urban_rural)
    results <- list(
      site = site_name,
      feeding_rate = feeding_rate,
      control = control,
      atsb_noirs = atsb_noirs,
      all_bells_whistles = atsb_bells_whistles
    )
    country_results[[site_name]] <- results
    print(j/length(country_site_files))
  }
  return(country_results)
}

par(las = 1)
plot(
  country_results$Bugiri_rural$all_bells_whistles$timestep/365 + 2000,
  country_results$Bugiri_rural$all_bells_whistles$n_detect_365_36499/
    country_results$Bugiri_rural$all_bells_whistles$n_365_36499,
  type = "l",
  lwd = 2,
  frame.plot = FALSE,
  col = "darkorchid4",
  xlab = "Year",
  ylab = "Prevalence",
  xlim = c(2020, 2025),
  ylim = c(0, 1)
)
grid()

lines(
  country_results$Bugiri_rural$atsb_noirs$timestep/365 + 2000,
  country_results$Bugiri_rural$atsb_noirs$n_detect_365_36499/
    country_results$Bugiri_rural$atsb_noirs$n_365_36499,
  lwd = 2,
  col = "mediumseagreen"
)

lines(
  country_results$Bugiri_rural$control$timestep/365 + 2000,
  country_results$Bugiri_rural$control$n_detect_365_36499/
    country_results$Bugiri_rural$control$n_365_36499,
  lwd = 2,
  col = "black"
)

legend("topleft", legend = c("Control", "ATSB replace IRS", "ATSB + IRS + IG2"), 
       col = c("black", "mediumseagreen", "darkorchid4"), lty = 1, bty = "n", 
       lwd = 2)

set.seed(1844)
params <- site_parameters(
  interventions = mali_sites[[4]]$interventions,
  demography = mali_sites[[4]]$demography,
  vectors = mali_sites[[4]]$vectors,
  seasonality = mali_sites[[4]]$seasonality,
  eir = mali_sites[[4]]$eir$eir[1],
  overrides = list(human_population = 5000)
)
params$bednet_timesteps[24] <- params$bednet_timesteps[23]+365
params$bednet_coverages[24] <- params$bednet_coverages[23]
params$bednet_dn0 <- rbind(params$bednet_dn0, params$bednet_dn0[23,])
params$bednet_rn <- rbind(params$bednet_rn, params$bednet_rn[23,])
params$bednet_rnm <- rbind(params$bednet_rnm, params$bednet_rnm[23,])
params$bednet_gamman[24] <- params$bednet_gamman[23]

params$spraying_timesteps[24] <- params$spraying_timesteps[23]
params$spraying_coverages[24] <- params$spraying_coverages[23]
params$spraying_ls_theta <- rbind(params$spraying_ls_theta, params$spraying_ls_theta[23,])
params$spraying_ls_gamma <- rbind(params$spraying_ls_gamma, params$spraying_ls_gamma[23,])
params$spraying_ks_theta <- rbind(params$spraying_ks_theta, params$spraying_ks_theta[23,])
params$spraying_ks_gamma <- rbind(params$spraying_ks_gamma, params$spraying_ks_gamma[23,])
params$spraying_ms_theta <- rbind(params$spraying_ms_theta, params$spraying_ms_theta[23,])
params$spraying_ms_gamma <- rbind(params$spraying_ms_gamma, params$spraying_ms_gamma[23,])

out2 <- run_simulation(timesteps = params$timesteps+5*365, parameters = params)

plot(out$timestep/365+2000,
     out$n_detect_730_3649/out$n_730_3649,
     type="l",
     lwd=2,
     frame.plot = F)
lines(out2$timestep/365+2000,
     out2$n_detect_730_3649/out$n_730_3649,
     lwd=2,
     col="darkorchid4")
