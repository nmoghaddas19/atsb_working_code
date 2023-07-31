library(foresite)
library(site)
library(malariasimulation)
library(grr)
library(cali)

pyr_nets <- read.csv("~/Documents/GitHub/pyrethroid_only_nets.csv")
ig2_nets <- read.csv("~/Documents/GitHub/pyrethroid_pyrrole_nets.csv")

uganda_districts <- data.frame(
  district = c("Bugiri", "Tororo", "Namutumba", "Butaleja", "Kibuku", "Budaka", 
               "Butebo", "Pallisa", "Serere", "Amolatar", "Kaberamaido", 
               "Kalaki", "Dokolo", "Alebtong", "Otuke", "Lira"),
  admin_1 = c("Bugiri", "Tororo", "Iganga", "Tororo", "Pallisa", "Pallisa",
              "Pallisa", "Pallisa", "Soroti", "Lira", "Kaberamaido", 
              "Kaberamaido", "Lira", "Lira", "Lira", "Lira"))

foo <- grr::matches(unique(uganda_districts$admin_1), UGA$sites$name_1, all.y = FALSE)
foo <- foo[order(foo$x),]
# 
# data <- list()
# for (i in 1:length(countries)) {
#   country_site_files <- countries[[names(countries)[i]]]
#   data[[names(countries)[i]]] <- run_country(country_site_files)
# }
country_results <- run_country_irs(UGA)

run_country_irs <- function(country_site_files) {
  country_name <- country_site_files$country
  country_results <- list()
  for (j in 1:length(foo$y)) {
    site <- single_site(UGA, foo$y[j])
    params <- site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir = site$eir$eir[1],
      overrides = list(human_population = 10000)
    )
    params$bednet_coverages[17:23] <- 0
    
    params$spraying_timesteps[24:26] <- c(params$spraying_timesteps[23] + 365,
                                          params$spraying_timesteps[23] + 365*2,
                                          params$spraying_timesteps[23] + 365*3)
    params$spraying_coverages[24:26] <- 0.6
    params$spraying_ls_theta <- rbind(params$spraying_ls_theta, params$spraying_ls_theta[23,],
                                       params$spraying_ls_theta[23,], params$spraying_ls_theta[23,])
    params$spraying_ls_gamma <- rbind(params$spraying_ls_gamma, params$spraying_ls_gamma[23,],
                                       params$spraying_ls_gamma[23,], params$spraying_ls_gamma[23,])
    params$spraying_ks_theta <- rbind(params$spraying_ks_theta, params$spraying_ks_theta[23,],
                                       params$spraying_ks_theta[23,], params$spraying_ks_theta[23,])
    params$spraying_ks_gamma <- rbind(params$spraying_ks_gamma, params$spraying_ks_gamma[23,],
                                       params$spraying_ks_gamma[23,], params$spraying_ks_gamma[23,])
    params$spraying_ms_theta <- rbind(params$spraying_ms_theta, params$spraying_ms_theta[23,],
                                       params$spraying_ms_theta[23,], params$spraying_ms_theta[23,])
    params$spraying_ms_gamma <- rbind(params$spraying_ms_gamma, params$spraying_ms_gamma[23,],
                                       params$spraying_ms_gamma[23,], params$spraying_ms_gamma[23,])
    params$spraying_coverages[17:23] <- 0.6
    
    
    EIR <- calibrate(
      parameters = params,
      target = site$prevalence$pfpr[20],
      summary_function = prev_at_baseline,
      tolerance = 0.02, 
      low = 100,
      high = 750,
      maxiter = 10
    )
    params <- set_equilibrium(parameters = params,
                              init_EIR = EIR)
    set.seed(j)
    control <- run_simulation(timesteps = params$timesteps + 5*365, 
                              parameters = params)
    
    atsb_stopirs_0.05 <- atsb_stopirs(0.05, site, EIR, j)
    atsb_stopirs_0.20 <- atsb_stopirs(0.10, site, EIR, j)
    atsb_stopirs_0.35 <- atsb_stopirs(0.35, site, EIR, j)
    
    bells_whistles_0.05 <- bells_whistles(0.05, site, EIR, j)
    bells_whistles_0.20 <- bells_whistles(0.10, site, EIR, j)
    bells_whistles_0.35 <- bells_whistles(0.35, site, EIR, j)

    site_name <- paste0(site$sites$name_1, "_", site$sites$urban_rural)
    results <- list(
      site = site_name,
      EIR = EIR,
      control = control,
      atsb_stopirs_0.05 = atsb_stopirs_0.05,
      atsb_stopirs_0.20 = atsb_stopirs_0.20,
      atsb_stopirs_0.35 = atsb_stopirs_0.35,
      bells_whistles_0.05 = bells_whistles_0.05,
      bells_whistles_0.20 = bells_whistles_0.20,
      bells_whistles_0.35 = bells_whistles_0.35
    )
    country_results[[site_name]] <- results
    print(j/length(foo$y))
  }
  return(country_results)
}
atsb_stopirs <- function(feeding_rate, site, EIR, j) {
  params_atsb_noirs <- site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    eir = site$eir$eir[1],
    overrides = list(human_population = 10000,
                     atsb = TRUE,
                     mu_atsb = c(feeding_rate, feeding_rate, feeding_rate))
  )
  params_atsb_noirs <- set_atsb(
    parameters = params_atsb_noirs,
    timesteps = (params_atsb_noirs$spraying_timesteps[23] + 365):(params_atsb_noirs$spraying_timesteps[23] + 4*365), 
    coverages = rep(1,1096)
  )
  res <- 71
  params_atsb_noirs$bednet_timesteps[24] <- params_atsb_noirs$spraying_timesteps[23] + 365
  params_atsb_noirs$bednet_coverages[24] <- max(params_atsb_noirs$bednet_coverages[21:23])
  params_atsb_noirs$bednet_dn0 <- rbind(params_atsb_noirs$bednet_dn0, rep(pyr_nets$dn0_med[res],3))
  params_atsb_noirs$bednet_rn <- rbind(params_atsb_noirs$bednet_rn, rep(pyr_nets$rn0_med[res],3))
  params_atsb_noirs$bednet_rnm <- rbind(params_atsb_noirs$bednet_rnm, params_atsb_noirs$bednet_rnm[23,])
  params_atsb_noirs$bednet_gamman[24] <- pyr_nets$gamman_med[res]
  params_atsb_noirs$bednet_coverages[17:23] <- 0
  
  params_atsb_noirs$spraying_coverages[17:23] <- 0.6

  params_atsb_noirs <- set_equilibrium(parameters = params_atsb_noirs,
                                       init_EIR = EIR)
  set.seed(j)
  atsb_noirs <- run_simulation(timesteps = params_atsb_noirs$timesteps + 5*365, 
                               parameters = params_atsb_noirs)
  return(atsb_noirs)
}
bells_whistles <- function(feeding_rate, site, EIR, j) {
  params_bells_whistles <- site_parameters(
    interventions = site$interventions,
    demography = site$demography,
    vectors = site$vectors,
    seasonality = site$seasonality,
    eir = site$eir$eir[1],
    overrides = list(human_population = 10000,
                     atsb = TRUE,
                     mu_atsb = c(feeding_rate, feeding_rate, feeding_rate))
  )
  params_bells_whistles <- set_atsb(
    parameters = params_bells_whistles,
    timesteps = (params_bells_whistles$spraying_timesteps[23] + 365):(params_bells_whistles$spraying_timesteps[23] + 4*365), 
    coverages = rep(1,1096)
  )
  
  res <- 71
  params_bells_whistles$bednet_timesteps[24] <- params_bells_whistles$spraying_timesteps[23] + 365
  params_bells_whistles$bednet_coverages[24] <- max(params_bells_whistles$bednet_coverages[21:23])
  params_bells_whistles$bednet_dn0 <- rbind(params_bells_whistles$bednet_dn0, rep(ig2_nets$dn0_med[res],3))
  params_bells_whistles$bednet_rn <- rbind(params_bells_whistles$bednet_rn, rep(ig2_nets$rn0_med[res],3))
  params_bells_whistles$bednet_rnm <- rbind(params_bells_whistles$bednet_rnm, params_bells_whistles$bednet_rnm[23,])
  params_bells_whistles$bednet_gamman[24] <- ig2_nets$gamman_med[res]*365
  params_bells_whistles$bednet_coverages[17:23] <- 0
  
  params_bells_whistles$spraying_timesteps[24:26] <- c(params_bells_whistles$spraying_timesteps[23] + 365,
                                                       params_bells_whistles$spraying_timesteps[23] + 365*2,
                                                       params_bells_whistles$spraying_timesteps[23] + 365*3)
  params_bells_whistles$spraying_coverages[24:26] <- 0.6
  params_bells_whistles$spraying_ls_theta <- rbind(params_bells_whistles$spraying_ls_theta, params_bells_whistles$spraying_ls_theta[23,],
                                                   params_bells_whistles$spraying_ls_theta[23,], params_bells_whistles$spraying_ls_theta[23,])
  params_bells_whistles$spraying_ls_gamma <- rbind(params_bells_whistles$spraying_ls_gamma, params_bells_whistles$spraying_ls_gamma[23,],
                                                   params_bells_whistles$spraying_ls_gamma[23,], params_bells_whistles$spraying_ls_gamma[23,])
  params_bells_whistles$spraying_ks_theta <- rbind(params_bells_whistles$spraying_ks_theta, params_bells_whistles$spraying_ks_theta[23,],
                                                   params_bells_whistles$spraying_ks_theta[23,], params_bells_whistles$spraying_ks_theta[23,])
  params_bells_whistles$spraying_ks_gamma <- rbind(params_bells_whistles$spraying_ks_gamma, params_bells_whistles$spraying_ks_gamma[23,],
                                                   params_bells_whistles$spraying_ks_gamma[23,], params_bells_whistles$spraying_ks_gamma[23,])
  params_bells_whistles$spraying_ms_theta <- rbind(params_bells_whistles$spraying_ms_theta, params_bells_whistles$spraying_ms_theta[23,],
                                                   params_bells_whistles$spraying_ms_theta[23,], params_bells_whistles$spraying_ms_theta[23,])
  params_bells_whistles$spraying_ms_gamma <- rbind(params_bells_whistles$spraying_ms_gamma, params_bells_whistles$spraying_ms_gamma[23,],
                                                   params_bells_whistles$spraying_ms_gamma[23,], params_bells_whistles$spraying_ms_gamma[23,])
  params_bells_whistles$spraying_coverages[17:23] <- 0.6
  
  params_bells_whistles <- set_equilibrium(parameters = params_bells_whistles,
                                           init_EIR = EIR)
  set.seed(j)
  atsb_bells_whistles <- run_simulation(timesteps = params_bells_whistles$timesteps + 5*365, 
                                        parameters = params_bells_whistles)
  return(atsb_bells_whistles)
}
saveRDS(object = country_results, 
        file = "~/Documents/uganda_IRS_3yearsATSB.RDS")
prev_at_baseline <- function(x) {
  baseline_timestep <- 19 * 365 + 182
  prev <- x[, "n_detect_730_3649"][baseline_timestep] / x[, "n_730_3649"][baseline_timestep]
  return(prev)
}
# 
# params$bednet_timesteps[24] <- params$bednet_timesteps[23] + 365
# params$bednet_coverages[24] <- max(params$bednet_coverages[21:23])
# params$bednet_dn0 <- rbind(params$bednet_dn0, params$bednet_dn0[23,])
# params$bednet_rn <- rbind(params$bednet_rn, params$bednet_rn[23,])
# params$bednet_rnm <- rbind(params$bednet_rnm, params$bednet_rnm[23,])
# params$bednet_gamman[24] <- params$bednet_gamman[23]

par(las = 1)
plot(
  country_results$Pallisa_rural$bells_whistles_0.20$timestep/365 + 2000,
  country_results$Pallisa_rural$bells_whistles_0.20$n_detect_730_3649/
    country_results$Pallisa_rural$bells_whistles_0.20$n_730_3649,
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
polygon(x = c(country_results$Pallisa_rural$bells_whistles_0.20$timestep/365 + 2000,
              rev(country_results$Pallisa_rural$bells_whistles_0.20$timestep/365 + 2000)),
        y = c(country_results$Pallisa_rural$bells_whistles_0.05$n_detect_730_3649/
                country_results$Pallisa_rural$bells_whistles_0.05$n_730_3649,
              rev(country_results$Pallisa_rural$bells_whistles_0.35$n_detect_730_3649/
                    country_results$Pallisa_rural$bells_whistles_0.35$n_730_3649)),
        border = FALSE,
        col = adjustcolor(col = "darkorchid4", alpha.f = 0.3))

lines(
  country_results$Pallisa_rural$atsb_stopirs_0.20$timestep/365 + 2000,
  country_results$Pallisa_rural$atsb_stopirs_0.20$n_detect_730_3649/
    country_results$Pallisa_rural$atsb_stopirs_0.20$n_730_3649,
  lwd = 2,
  col = "mediumseagreen"
)

polygon(x = c(country_results$Pallisa_rural$atsb_stopirs_0.20$timestep/365 + 2000,
              rev(country_results$Pallisa_rural$atsb_stopirs_0.20$timestep/365 + 2000)),
        y = c(country_results$Pallisa_rural$atsb_stopirs_0.05$n_detect_730_3649/
                country_results$Pallisa_rural$atsb_stopirs_0.05$n_730_3649,
              rev(country_results$Pallisa_rural$atsb_stopirs_0.35$n_detect_730_3649/
                    country_results$Pallisa_rural$atsb_stopirs_0.35$n_730_3649)),
        border = FALSE,
        col = adjustcolor(col = "mediumseagreen", alpha.f = 0.3))

lines(
  country_results$Pallisa_rural$control$timestep/365 + 2000,
  country_results$Pallisa_rural$control$n_detect_730_3649/
    country_results$Pallisa_rural$control$n_730_3649,
  lwd = 2,
  col = "black"
)

legend("topleft", legend = c("Control", "ATSB replace IRS", "ATSB + IRS + IG2"), 
       col = c("black", "mediumseagreen", "darkorchid4"), lty = 1, bty = "n", 
       lwd = 2)

par(las = 1)
plot(
  control$timestep/365 + 2000,
  control$n_detect_730_3649/
    control$n_730_3649,
  type = "l",
  lwd = 2,
  frame.plot = FALSE,
  col = "darkorchid4",
  xlab = "Year",
  ylab = "Prevalence",
  xlim = c(2010, 2025),
  ylim = c(0, 1)
)
grid()

set.seed(1844)
params <- site_parameters(
  interventions = mali_sites[[4]]$interventions,
  demography = mali_sites[[4]]$demography,
  vectors = mali_sites[[4]]$vectors,
  seasonality = mali_sites[[4]]$seasonality,
  eir = mali_sites[[4]]$eir$eir[1],
  overrides = list(human_population = 10000)
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
