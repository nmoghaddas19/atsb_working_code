ghana_data <- readRDS("~/Documents/ghana_IRS_3yearsATSB.RDS")
EIRs <- unlist(lapply(ghana_data, "[[", 2))

library(foresite)
library(site)
library(malariasimulation)
library(grr)
library(cali)

pyr_nets <- read.csv("~/Documents/GitHub/pyrethroid_only_nets.csv")
ig2_nets <- read.csv("~/Documents/GitHub/pyrethroid_pyrrole_nets.csv")

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

foo <- grr::matches(unique(ghana_districts$admin_1), GHA$sites$name_1, all.y = FALSE)
foo <- foo[order(foo$x),]
# 
# data <- list()
# for (i in 1:length(countries)) {
#   country_site_files <- countries[[names(countries)[i]]]
#   data[[names(countries)[i]]] <- run_country(country_site_files)
# }
country_results <- run_country_irs(GHA)

run_country_irs <- function(country_site_files) {
  country_name <- country_site_files$country
  country_results <- list()
  for (j in 1:length(foo$y)) {
    message(j/length(foo$y))
    site <- single_site(GHA, foo$y[j])
    site$interventions$rtss_cov[1:23] <- 0
    params <- site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir = site$eir$eir[1],
      overrides = list(human_population = 10000)
    )
    
    EIR <- EIRs[j]
    
    res <- 71
    params$bednet_timesteps[24] <- params$spraying_timesteps[23] + 365
    params$bednet_coverages[24] <- max(params$bednet_coverages[21:23])
    params$bednet_dn0 <- rbind(params$bednet_dn0, rep(pyr_nets$dn0_med[res],3))
    params$bednet_rn <- rbind(params$bednet_rn, rep(pyr_nets$rn0_med[res],3))
    params$bednet_rnm <- rbind(params$bednet_rnm, params$bednet_rnm[23,])
    params$bednet_gamman[24] <- pyr_nets$gamman_med[res]
    
    params$spraying_coverages[1:23] <- 0
    
    params <- set_equilibrium(parameters = params,
                              init_EIR = EIR)
    set.seed(j)
    control <- run_simulation(timesteps = params$timesteps + 5*365, 
                              parameters = params)
    
    atsb_stopirs_0.05 <- atsb(0.05, site, EIR, j, params) # this is really adding atsb to control, stopirs is
    atsb_stopirs_0.20 <- atsb(0.10, site, EIR, j, params) # wrong name but keeping it so code runs easily
    atsb_stopirs_0.35 <- atsb(0.35, site, EIR, j, params)
    
    bells_whistles_0.05 <- bells_whistles(0.05, site, EIR, j, params)
    bells_whistles_0.20 <- bells_whistles(0.10, site, EIR, j, params)
    bells_whistles_0.35 <- bells_whistles(0.35, site, EIR, j, params)
    
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
atsb <- function(feeding_rate, site, EIR, j, params) {
  params_atsb <- set_atsb(
    parameters = params,
    timesteps = (params$spraying_timesteps[23] + 365):(params$spraying_timesteps[23] + 4*365), 
    coverages = rep(1,1096)
  )
  
  params_atsb$atsb = TRUE
  params_atsb$mu_atsb = c(feeding_rate, feeding_rate, feeding_rate)
  
  params_atsb <- set_equilibrium(parameters = params_atsb,
                                 init_EIR = EIR)
  set.seed(j)
  atsb_noirs <- run_simulation(timesteps = params_atsb$timesteps + 5*365, 
                               parameters = params_atsb)
  return(atsb_noirs)
}
bells_whistles <- function(feeding_rate, site, EIR, j, params) {
  
  params_bells_whistles <- set_atsb(
    parameters = params,
    timesteps = (params$spraying_timesteps[23] + 365):(params$spraying_timesteps[23] + 4*365), 
    coverages = rep(1,1096)
  )
  
  params_bells_whistles$atsb = TRUE
  params_bells_whistles$mu_atsb = c(feeding_rate, feeding_rate, feeding_rate)
  
  res <- 71
  params_bells_whistles$bednet_timesteps[24] <- params_bells_whistles$spraying_timesteps[23] + 365
  params_bells_whistles$bednet_coverages[24] <- max(params_bells_whistles$bednet_coverages[21:23])
  params_bells_whistles$bednet_dn0[24,] <-  rep(ig2_nets$dn0_med[res],3)
  params_bells_whistles$bednet_rn[24,] <-  rep(ig2_nets$rn0_med[res],3)
  params_bells_whistles$bednet_rnm <- rbind(params_bells_whistles$bednet_rnm, params_bells_whistles$bednet_rnm[23,])
  params_bells_whistles$bednet_gamman[24] <- ig2_nets$gamman_med[res]*365
  
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
