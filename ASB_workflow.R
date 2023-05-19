################################################################################
# In this script, I have written a workflow that takes a dataframe of the      #
# sort we are expecting from the ATSB trials and generates predictions with    #
# confidence intervals based on malariasimulation. The first step is to        #
# generate estimates of dyed fraction and uncertainty from the data. The       #
# second step is to convert the dyed fraction to daily feeding rates with.     #
# uncertainty. The third step is to plug those bounds into malariasim to       #
# generate predictions.                                                        #
################################################################################

setwd("~/GitHub/")
# load in packages
library(malariasimulation)
library(ICDMM)
library(site)
library(foresite)
library(lubridate)
library(RColorBrewer)
library(lme4)
library(dplyr)
library(cali)
library(umbrella)

# functions
jensen.logit.adjust <-
  function(p, V, method = "mcculloch", inverse = FALSE) {
    if(method == "mcculloch" & inverse) {
      method <- "zeger"
      warning("The McCulloch method can't be used when inverse = TRUE. Changing to Zeger.")
    }
    stopifnot(!(method == "mcculloch" & inverse))
    Beta <- qlogis(p)
    if(method == "mcculloch") {
      return(plogis(Beta - 0.5 * V * tanh(Beta * (1 + 2 * exp(-0.5 * V))/6)))
    }
    if(method == "zeger") {
      if (inverse) {
        plogis(Beta * sqrt(256 * V / (75 * pi^2) + 1))
      } else {
        plogis(Beta/sqrt(1 + ((16 * sqrt(3))/(15 * pi))^2 * V))
      }
    }
  }

# load in data
dat.all <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")
# cleaning and checking ####
names(dat.all) <- tolower(names(dat.all))
dim(dat.all)
table(table(dat.all$form.sampleid.check))
dat.all <- dat.all[!duplicated(dat.all$form.sampleid.check), ]
table(table(dat.all$form.sampleid.check))
dat.all$mosquito_species[dat.all$an..gambiae %in% 1] <- "gambiae" # make species a variable
dat.all$mosquito_species[dat.all$an..funestus %in% 1] <- "funestus" # make species a variable
dat.all$mosquito_species <- factor(dat.all$mosquito_species) # make species a factor
dat.all$location <- factor(dat.all$indoor, 0:1, c("Outdoors", "Indoors")) # make location a factor
dat.all$anoph_sex <- factor(dat.all$male, 0:1, c("female", "male")) # make sex a factor
dat.all$n.asb <- dat.all$asbs.deployed 
table(dat.all$n.asb, exclude = NULL)
dat.all$date.id <- as.Date(dat.all$day.of.date.trap, "%d-%b-%y")
dat.all$month_caught <- factor(month(dat.all$date.id))
dat.all$collection_date <- factor(dat.all$date.id)
dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 2] <-
  paste0(dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 2], ".23")
dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 3] <-
  paste0(dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 3], ".23")
dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 3] <-
  paste0(dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 3], ".32")
dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 2] <-
  paste0(dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 2], ".32")

dat.all$an..gambiae <- dat.all$an..funestus <- dat.all$indoor <- dat.all$outdoor <-
  dat.all$male <- dat.all$female <- dat.all$asbs.deployed <-
  dat.all$day.of.date.trap <- dat.all$cluster.code <- NULL
sum(is.na(dat.all))
colSums(is.na(dat.all))
table(dat.all$form.bloodfed, exclude = NULL)
dat.all <-
  dat.all[dat.all$mosquito_species %in% c("gambiae", "funestus") &
            dat.all$anoph_sex %in% "female", ]
table(dat.all$cluster)
table(dat.all$mosquito_species, exclude = NULL)
table(dat.all$collection_method, exclude = NULL)
table(dat.all$anoph_sex, exclude = NULL)
range(dat.all$date.id)
tapply(dat.all$date.id, dat.all$phase, range)
table(dat.all$cluster, dat.all$n.asb)
table(dat.all$cluster, dat.all$n.asb, dat.all$phase)
dat.all <-
  dat.all[order(dat.all$n.asb, dat.all$cluster, dat.all$collection_date,
                dat.all$mosquito_species, dat.all$location), ]
dat.all$cluster <- factor(dat.all$cluster, unique(dat.all$cluster))
dat.all$n.asb <- factor(dat.all$n.asb, 2:3, c("2 stations", "3 stations"))
dat.all$hh <- factor(dat.all$form.collectionid.check)
dat.all$form.collectionid.check <- NULL
table(form.bloodfed = dat.all$form.bloodfed, form.dyefed = dat.all$form.dyefed, 
      exclude = NULL)
prop.table(table(form.bloodfed = dat.all$form.bloodfed, 
                 form.dyefed = dat.all$form.dyefed, exclude = NULL))
dat.all$positive <- dat.all$form.dyefed
dat.all$form.dyefed <- NULL
dat.all$cluster.day <- factor(paste(dat.all$cluster, dat.all$date.id, sep = ":"))

# stop ####

# fit regression model for each species 
sp <- c("funestus")
dat <- droplevels(dat.all[dat.all$mosquito_species  %in% sp, ])
term <- "cluster"
rest.of.formula <- "positive ~ (1 | collection_date) + (1|hh)"
form <- paste(rest.of.formula, term, sep = " + ")                             
fit <- glmer(form, family = "binomial", data = dat, 
             control = glmerControl(optimizer = "bobyqa"))
print(summary(fit))

fit0 <- update(fit, ~ . -1)
est.tab.bias <-
  cbind(fixef(fit0),
        confint(fit0, method = "Wald")[paste0(term, levs), ])
# adjust predictions for Jensen's inequality
est.tab_fun <-
  jensen.logit.adjust(p = plogis(est.tab.bias), V = sum(unlist(VarCorr(fit0))))

# dat.all$form.bloodfed <- factor(dat.a;;)
sp <- c("gambiae")
dat <- droplevels(dat.all[dat.all$mosquito_species  %in% sp, ])
term <- "cluster"
rest.of.formula <- "positive ~ (1 | collection_date)"
form <- paste(rest.of.formula, term, sep = " + ")                             
fit <- glmer(form, family = "binomial", data = dat, 
             control = glmerControl(optimizer = "bobyqa"))
print(summary(fit))

fit0 <- update(fit, ~ . -1)
est.tab.bias <-
  cbind(fixef(fit0),
        confint(fit0, method = "Wald")[paste0(term, levs), ])
# adjust predictions for Jensen's inequality
est.tab_gamb <-
  jensen.logit.adjust(p = plogis(est.tab.bias), V = sum(unlist(VarCorr(fit0))))

dyed_fraction_bounds <- cbind(est.tab_fun, est.tab_gamb)
colnames(dyed_fraction_bounds)[c(1,4)] <- c("fun", "gamb")

# convert to daily feeding rate
dyed_fraction <- read.csv(file="~/Documents/GitHub/atsb_working_code/dyed_fraction.csv")[1:200,]
dyed_fraction <- apply(dyed_fraction, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975))[,7:31]
upper_fun <- approx(x = dyed_fraction[1,], y = seq(1,50,2)/100,
                xout = dyed_fraction_bounds[,3])$y
lower_fun <- approx(x = dyed_fraction[2,], y = seq(1,50,2)/100, 
                xout = dyed_fraction_bounds[,2])$y

upper_gamb <- approx(x = dyed_fraction[1,], y = seq(1,50,2)/100,
                     xout = dyed_fraction_bounds[,6])$y
lower_gamb <- approx(x = c(0,dyed_fraction[2,]), y = c(0,seq(1,50,2)/100), 
                     xout = dyed_fraction_bounds[,5])$y

feed_rates <- cbind(control = numeric(10), lower_fun, upper_fun, lower_gamb, upper_gamb)
rownames(feed_rates) <- names(fixef(fit0))

# forecast predictions using malariasimulation
# !!! TWO APPROCAHES !!! 
# OPTION ONE: using historical assumptions in foresite site files ####
zambia <- ZMB
zambia$sites$name_1
western_rural <- single_site(zambia, 18)
# plot(western_rural$prevalence$year, western_rural$prevalence$pfpr, type="l", lwd=2, 
#      frame.plot = F, ylim = c(0,1))
# western_rural$vectors

out_data <- list()
for (j in 1:2) {
  bound <- c("lower", "upper")[j]
  for (i in 1:nrow(feed_rates)) {
    cluster <- names(fixef(fit0))[i]
    western_rural_params <- site_parameters(
      interventions = western_rural$interventions,
      demography = western_rural$demography,
      vectors = western_rural$vectors,
      seasonality = western_rural$seasonality,
      eir = western_rural$eir$eir[1],
      overrides = list(human_population = 10000,
                       mu_atsb = c(feed_rates[i,j+3], feed_rates[i,j+1], feed_rates[i,j+3]))
    )
    western_rural_params <- set_atsb(parameters = western_rural_params,
                                     timesteps = (21*365+60):(21*365+148), 
                                     coverages = rep(1,148-59))
    name <- paste(cluster, bound, sep = "_")
    out_data[[name]] <- run_simulation(timesteps = western_rural_params$timesteps,
                                          parameters = western_rural_params)
    print(i/nrow(feed_rates))
  }
  print(j)
}
western_rural_params <- site_parameters(
  interventions = western_rural$interventions,
  demography = western_rural$demography,
  vectors = western_rural$vectors,
  seasonality = western_rural$seasonality,
  eir = western_rural$eir$eir[1],
  overrides = list(human_population = 10000,
                   mu_atsb = c(0,0,0))
)
western_rural_params <- set_atsb(parameters = western_rural_params,
                                 timesteps = (21*365+60):(21*365+148), 
                                 coverages = rep(1,148-59))
for (i in 1:10) {
  if (i==1) {
    out_data[["Control"]] <- run_simulation(timesteps = western_rural_params$timesteps,
                                       parameters = western_rural_params)
  } else {
    name <- paste0("Control","_ ", i)
    out_data[[name]] <- run_simulation(timesteps = western_rural_params$timesteps,
                                       parameters = western_rural_params)
  }
  print(i/10)
}

plot(out_data$Control$timestep/365 +2000, 
     out_data$Control$total_M_arabiensis+out_data$Control$total_M_funestus+out_data$Control$total_M_gambiae,
     type="l", frame.plot = F, xlim = c(2019,2022),
     xlab = "Year", ylab = "Total mosquitoes")
# i=2
# total_M_lower <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
# total_M_upper <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
# polygon(c(out_data[[i]]$timestep/365 +2000, rev(out_data[[i+10]]$timestep/365 +2000)),
#         c(total_M_lower,rev(total_M_upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
#         border = F)
for (i in 1:10) {
  total_M_lower <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
  total_M_upper <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
  total_M_control <- out_data[[i+20]]$total_M_arabiensis + out_data[[i+20]]$total_M_funestus + out_data[[i+20]]$total_M_gambiae
  lines(out_data[[i]]$timestep/365+2000, total_M_lower, col = brewer.pal(10, "RdYlBu")[2])
  lines(out_data[[i+10]]$timestep/365+2000, total_M_upper, col = brewer.pal(10, "RdYlBu")[9])
  lines(out_data[[i+20]]$timestep/365+2000, total_M_control, col="black")
}
lower <- matrix(0,nrow = 8395, ncol=10)
upper <- matrix(0,nrow = 8395, ncol=10)
control <- matrix(0,nrow = 8395, ncol=10)
for (i in 1:10) {
  lower[,i] <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
  upper[,i] <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
  control[,i] <- out_data[[i+20]]$total_M_arabiensis + out_data[[i+20]]$total_M_funestus + out_data[[i+20]]$total_M_gambiae
}
lower <- apply(lower, MARGIN = 1, FUN = mean)
upper <- apply(upper, MARGIN = 1, FUN = mean)
control <- apply(control, MARGIN = 1, FUN = mean)
plot(out_data$Control$timestep/365 +2000, 
     control,
     type="l", frame.plot = F, lwd=2, xlim = c(2019,2022),
     xlab = "Year", ylab = "Total mosquitoes", ylim=c(0,3000000))
polygon(c(out_data$Control$timestep/365 +2000, rev(out_data$Control$timestep/365 +2000)),
        c(lower,rev(upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
        border = F)

plot(out_data$Control$timestep/365 +2000, 
     out_data$Control$n_detect_730_3649/out_data$Control$n_730_3649,
     type="l", frame.plot = F, xlim = c(2019,2022), ylim=c(0,0.8),
     xlab = "Year", ylab = "PfPr2-10")
# i=2
# total_M_lower <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
# total_M_upper <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
# polygon(c(out_data[[i]]$timestep/365 +2000, rev(out_data[[i+10]]$timestep/365 +2000)),
#         c(total_M_lower,rev(total_M_upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
#         border = F)
for (i in 1:10) {
  pfpr_lower <- out_data[[i]]$n_detect_730_3649/out_data[[i]]$n_730_3649
  pfpr_upper <- out_data[[i+10]]$n_detect_730_3649/out_data[[i+10]]$n_730_3649
  pfpr_control <- out_data[[i+20]]$n_detect_730_3649/out_data[[i+20]]$n_730_3649
  lines(out_data[[i]]$timestep/365+2000, pfpr_lower, col = brewer.pal(10, "RdYlBu")[2])
  lines(out_data[[i+10]]$timestep/365+2000, pfpr_upper, col = brewer.pal(10, "RdYlBu")[9])
  lines(out_data[[i+20]]$timestep/365+2000, pfpr_control, col="black")
}
lower <- matrix(0,nrow = 8395, ncol=10)
upper <- matrix(0,nrow = 8395, ncol=10)
control <- matrix(0,nrow = 8395, ncol=10)
for (i in 1:10) {
  lower[,i] <- out_data[[i]]$n_detect_730_3649/out_data[[i]]$n_730_3649
  upper[,i] <- out_data[[i+10]]$n_detect_730_3649/out_data[[i+10]]$n_730_3649
  control[,i] <- out_data[[i+20]]$n_detect_730_3649/out_data[[i+20]]$n_730_3649
}
lower <- apply(lower, MARGIN = 1, FUN = mean)
upper <- apply(upper, MARGIN = 1, FUN = mean)
control <- apply(control, MARGIN = 1, FUN = mean)
plot(out_data$Control$timestep/365 +2000, 
     control,
     type="l", frame.plot = F, lwd=2, xlim = c(2019,2022),
     xlab = "Year", ylab = "PfPr2-10", ylim=c(0,0.2))
polygon(c(out_data$Control$timestep/365 +2000, rev(out_data$Control$timestep/365 +2000)),
        c(lower,rev(upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
        border = F)

# stop ####

# OPTION TWO: calibrating to baseline prevalence and predicting forwards ####
zambia <- ZMB
western_rural <- single_site(zambia, 18)

# int <- western_rural$interventions
# int <- rbind(int, int[rep(23,3),])
# int$year[24:26] <- c(2023:2025)
# int <- int[c(18,21,24),]
# 
# demog <- western_rural$demography
# demog |>
#   filter(year>=2017) -> demog
# vectors <- western_rural$vectors
# western_rural$seasonality
# western_rural$eir
# 
# western_rural_params <- site_parameters(
#   interventions = int,
#   demography = demog,
#   vectors = western_rural$vectors,
#   seasonality = western_rural$seasonality,
#   eir = 340,
#   overrides = list(human_population = 10000,
#                    mu_atsb = c(0,0,0))
# )

# starting simulation in 2014, ITN distributions in 2015 and 2018, baseline at 2021
data <- data_frame(Cluster = 1:10, 
                   feed_rate_arab = rnorm(10, mean=0.14, sd=0.05),
                   feed_rate_arab_lo = rnorm(10, mean=0.08, sd=0.04),
                   feed_rate_arab_hi = rnorm(10, mean=0.21, sd=0.04),
                   feed_rate_fun = rnorm(10, mean=0.14, sd=0.05),
                   feed_rate_fun_lo = rnorm(10, mean=0.08, sd=0.04),
                   feed_rate_fun_hi = rnorm(10, mean=0.21, sd=0.04),
                   feed_rate_gamb = rnorm(10, mean=0.14, sd=0.05),
                   feed_rate_gamb_lo = rnorm(10, mean=0.08, sd=0.04),
                   feed_rate_gamb_hi = rnorm(10, mean=0.21, sd=0.04),
                   pfpr_baseline = rnorm(10, mean=0.2, sd=0.05),
                   baseline_time = round(rnorm(10, mean=7*365+50, sd=20), digits = 0),
                   itn_timestep_1 = rep(365-30, 10),
                   itn_coverage_1 = rep(western_rural$interventions$itn_use[16], 10),
                   resistance_1 = rep(western_rural$interventions$pyrethroid_resistance[16], 10),
                   itn_timestep_2 = rep(4*365-30, 10),
                   itn_coverage_2 = rep(western_rural$interventions$itn_use[19], 10),
                   resistance_2 = rep(western_rural$interventions$pyrethroid_resistance[19], 10),
                   itn_timestep_3 = rep(7*365-30, 10),
                   itn_coverage_3 = rep(western_rural$interventions$itn_use[22], 10),
                   resistance_3 = rep(0.52, 10),
                   retention = rep(5*365, 10),
                   irs_timestep_1 = rep(365-30, 10),
                   irs_coverage_1 = rep(western_rural$interventions$irs_cov[16], 10),
                   irs_timestep_2 = rep(4*365-30, 10),
                   irs_coverage_2 = rep(western_rural$interventions$irs_cov[19], 10),
                   irs_timestep_3 = rep(7*365-30, 10),
                   irs_coverage_3 = rep(western_rural$interventions$irs_cov[22], 10), 
                   ft_1 = rep(western_rural$interventions$tx_cov[16], 10),
                   prop_act_1 = rep(western_rural$interventions$prop_act[16], 10),
                   ft_2 = rep(western_rural$interventions$tx_cov[19], 10),
                   prop_act_2 = rep(western_rural$interventions$prop_act[19], 10),
                   ft_3 = rep(western_rural$interventions$tx_cov[22], 10),
                   prop_act_3 = rep(western_rural$interventions$prop_act[22], 10),
                   arab = rep(western_rural$vectors$prop[1], 10),
                   fun = rep(western_rural$vectors$prop[2], 10),
                   gamb = rep(western_rural$vectors$prop[3], 10),
                   g0 = rep(western_rural$seasonality$g0, 10),
                   g1 = rep(western_rural$seasonality$g1, 10),
                   g2 = rep(western_rural$seasonality$g2, 10),
                   g3 = rep(western_rural$seasonality$g3, 10),
                   h1 = rep(western_rural$seasonality$h1, 10),
                   h2 = rep(western_rural$seasonality$h2, 10),
                   h3 = rep(western_rural$seasonality$h3, 10)
                   
                   )

cluster_params <- get_parameters(list(
  human_population = 10000,
  model_seasonality = TRUE,
  g0 = data$g0[i],
  g = c(data$g1[i], data$g2[i], data$g3[i]),
  h = c(data$h1[i], data$h2[i], data$h3[i]),
  individual_mosquitoes = FALSE
))

cluster_params <- set_species(parameters = cluster_params,
                              species = list(arab_params, fun_params, gamb_params),
                              proportions = c(data$arab[i], data$fun[i], data$gamb[i]))

cluster_params <- set_drugs(parameters = cluster_params, drugs = list(AL_params))
cluster_params <- set_clinical_treatment(parameters = cluster_params,
                                         drug = 1,
                                         timesteps = c(0,3,6)*365,
                                         coverages = c(data$ft_1[i]*data$prop_act_1[i],
                                                       data$ft_2[i]*data$prop_act_2[i],
                                                       data$ft_3[i]*data$prop_act_3[i]))
pyrethroid_params <- read.csv("pyrethroid_only_nets.csv")
dn0_1 <- pyrethroid_params$dn0_med[data$resistance_1[i]*100+1]
dn0_2 <- pyrethroid_params$dn0_med[data$resistance_2[i]*100+1]
dn0_3 <- pyrethroid_params$dn0_med[data$resistance_3[i]*100+1]
rn0_1 <- pyrethroid_params$rn0_med[data$resistance_1[i]*100+1]
rn0_2 <- pyrethroid_params$rn0_med[data$resistance_2[i]*100+1]
rn0_3 <- pyrethroid_params$rn0_med[data$resistance_3[i]*100+1]
gamman_1 <- pyrethroid_params$gamman_med[data$resistance_1[i]*100+1]
gamman_2 <- pyrethroid_params$gamman_med[data$resistance_2[i]*100+1]
gamman_3 <- pyrethroid_params$gamman_med[data$resistance_3[i]*100+1]

cluster_params <- set_bednets(
  parameters = cluster_params,
  timesteps = c(data$itn_timestep_1[i], data$itn_timestep_2[i], data$itn_timestep_3[i]),
  coverages = c(data$itn_coverage_1[i], data$itn_coverage_2[i], data$itn_coverage_3[i]),
  retention = data$retention[i],
  dn0 = matrix(c(dn0_1, dn0_2, dn0_3,
                 dn0_1, dn0_2, dn0_3,
                 dn0_1, dn0_2, dn0_3), nrow =3, ncol=3),
  rn = matrix(c(rn0_1, rn0_2, rn0_3,
                rn0_1, rn0_2, rn0_3,
                rn0_1, rn0_2, rn0_3), nrow =3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow = 3, ncol = 3),
  gamman = c(gamman_1, gamman_2, gamman_3)*365
)
ls_theta_1 = western_rural$interventions$ls_theta[16]
ls_gamma_1 = western_rural$interventions$ls_gamma[16]
ks_theta_1 = western_rural$interventions$ks_theta[16]
ks_gamma_1 = western_rural$interventions$ks_gamma[16]
ms_theta_1 = western_rural$interventions$ms_theta[16]
ms_gamma_1 = western_rural$interventions$ms_gamma[16]
ls_theta_2 = western_rural$interventions$ls_theta[19]
ls_gamma_2 = western_rural$interventions$ls_gamma[19]
ks_theta_2 = western_rural$interventions$ks_theta[19]
ks_gamma_2 = western_rural$interventions$ks_gamma[19]
ms_theta_2 = western_rural$interventions$ms_theta[19]
ms_gamma_2 = western_rural$interventions$ms_gamma[19]
ls_theta_3 = western_rural$interventions$ls_theta[22]
ls_gamma_3 = western_rural$interventions$ls_gamma[22]
ks_theta_3 = western_rural$interventions$ks_theta[22]
ks_gamma_3 = western_rural$interventions$ks_gamma[22]
ms_theta_3 = western_rural$interventions$ms_theta[22]
ms_gamma_3 = western_rural$interventions$ms_gamma[22]

cluster_params <- set_spraying(
  parameters = cluster_params,
  timesteps = c(data$irs_timestep_1[i], data$irs_timestep_2[i], data$irs_timestep_3[i]),
  coverages = c(data$irs_coverage_1[i], data$irs_coverage_2[i], data$irs_coverage_3[i]),
  ls_theta = matrix(c(ls_theta_1, ls_theta_2, ls_theta_3,
                      ls_theta_1, ls_theta_2, ls_theta_3,
                      ls_theta_1, ls_theta_2, ls_theta_3), nrow =3, ncol=3),
  ls_gamma = matrix(c(ls_gamma_1, ls_gamma_2, ls_gamma_3,
                      ls_gamma_1, ls_gamma_2, ls_gamma_3,
                      ls_gamma_1, ls_gamma_2, ls_gamma_3), nrow =3, ncol=3),
  ks_theta = matrix(c(ks_theta_1, ks_theta_2, ks_theta_3,
                      ks_theta_1, ks_theta_2, ks_theta_3,
                      ks_theta_1, ks_theta_2, ks_theta_3), nrow =3, ncol=3),
  ks_gamma = matrix(c(ks_gamma_1, ks_gamma_2, ks_gamma_3,
                      ks_gamma_1, ks_gamma_2, ks_gamma_3,
                      ks_gamma_1, ks_gamma_2, ks_gamma_3), nrow =3, ncol=3),
  ms_theta = matrix(c(ms_theta_1, ms_theta_2, ms_theta_3,
                      ms_theta_1, ms_theta_2, ms_theta_3,
                      ms_theta_1, ms_theta_2, ms_theta_3), nrow =3, ncol=3),
  ms_gamma = matrix(c(ms_gamma_1, ms_gamma_2, ms_gamma_3,
                      ms_gamma_1, ms_gamma_2, ms_gamma_3,
                      ms_gamma_1, ms_gamma_2, ms_gamma_3), nrow =3, ncol=3)
)

cluster_params$timesteps <- 10*365
prev_at_baseline <- function(x) {
  baseline_timestep <- data$baseline_time[i]
  prev <- x[, "n_detect_730_3650"][baseline_timestep] / x[, "n_730_3650"][baseline_timestep]
  return(prev)
}

set.seed(123)
EIR <- calibrate(
  parameters = cluster_params,
  target = data$pfpr_baseline[i],
  summary_function = prev_at_baseline, #using own summary function at required baseline survey timestep
  tolerance = 0.01,
  low = 0.001,
  high = 350
)
cluster_params <- set_equilibrium(
  parameters = cluster_params,
  init_EIR = EIR
)
out <- run_simulation(timesteps = 10*365,
                      parameters = cluster_params)

par(las=1)
plot(out$timestep/365+2014, out$n_detect_730_3650/out$n_730_3650,
     type = "l", lwd=2, xlab="Year", ylab="PfPr2-10", frame.plot = F,
     ylim=c(0,1), xlim=c(2014,2024))
grid()
abline(v=itn_distributions/365+2015, lty=2)
points(x= data$baseline_time[i]/365+2014, y = data$pfpr_baseline[i],
       cex=3.5, col="mediumseagreen", pch=20)
# stop ####




