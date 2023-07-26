################################################################################
## The purpose of this script is to use learn maximum likelihood estimation   ##
## by using using it to fit the feeding rate parameter to the observed count  ##
## data from Mali in 2017. A similar exercise was carried out in Fraser et al.##
## 2020 using a different model and a different method. Hopefully this will   ##
## also give us an idea of how differently they are estimating catch numbers. ##
################################################################################

library(malariasimulation)
library(dplyr)
library(foresite)
library(site)

# load in the count data
cdc_2017 <- read.csv("~/Documents/GitHub/atsb_working_code/DB CDC Malaise PSC catches 2017/CDC-Table 1.csv")

# load control malariasim run
malariasim_control <- readRDS("~/Documents/GitHub/atsb_working_code/out_mali.RDS")$Kayes_rural_data

# find scaling factor 
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=2.44*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> cdc_2017_con
scaler <- max(cdc_2017_con$Mean)/max(malariasim_control$total_M_gambiae[(17*365):(18*365)])

# visualise
plot(malariasim_control$timestep/365+2000, malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018))
lines(2017+(cdc_2017_con$Month)/12, cdc_2017_con$Mean, lwd=2, col=1, lty=2)

# calculate sum of least squares
sum_of_squares <- c()

for (i in 1:80) {
  true_values <- cdc_2017_con$Mean 
  y <- malariasim_control$total_M_gambiae*scaler
  x <- malariasim_control$timestep/365 + 2000 + i/365
  model_values <- approx(x, y, xout= 2017 + c(4:12)/12)
  sum_of_squares[i] <- sum((true_values - model_values$y)^2)
}
plot(1:80, sum_of_squares, type="l", lwd=2, frame.plot = F, xlab="Days shifted")

days_shifted <- which(sum_of_squares == min(sum_of_squares))

# lets see the data with that shift
par(las=1,  = )
plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018),
     xlab="Year", ylab="Population", cex.lab=1.2)
polygon(c(out_bounds[[2]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[3]]$timestep/365+2000+days_shifted/365)),
        c(out_bounds[[2]]$total_M_gambiae*scaler, rev(out_bounds[[3]]$total_M_gambiae)*scaler),
        col = adjustcolor("dodgerblue", alpha.f = 0.5), border = FALSE)
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> cdc_2017_con
lines(2017+(cdc_2017_con$Month)/12, cdc_2017_con$Mean, lwd=2, col=1, lty=2)
arrows(2017+(cdc_2017_con$Month)/12,
       cdc_2017_con$Mean-cdc_2017_con$CI,
       2017+(cdc_2017_con$Month)/12,
       cdc_2017_con$Mean+cdc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.05,
       col=1)
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> cdc_2017_exp
lines(2017+(cdc_2017_exp$Month)/12, cdc_2017_exp$Mean, lwd=2, col=2, lty=2)
arrows(2017+(cdc_2017_exp$Month)/12, 
       cdc_2017_exp$Mean-cdc_2017_exp$CI,
       2017+(cdc_2017_exp$Month)/12,
       cdc_2017_exp$Mean+cdc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=2)
legend(x="topleft", legend=c("Control data", "ATSB data", "Control model", "ATSB model"), 
       col=c(1,2,1,"dodgerblue"), lwd=2, lty=c(2,2,1,1), bty="n")
title("CDC traps")

# bootstrap CIs ----

d_con <- matrix(0, nrow=7, ncol=3)
for (i in 6:12) {
  cdc_2017 |>
    filter(Month == i & Experimental.or.control == "Con.") -> t
  x <- matrix(rep(t$tot.f, 5000), byrow = T, ncol = length(t$tot.f))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,7,T))-mean(x)})
  d_con[i-5,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$tot.f)
}

d_exp <- matrix(0, nrow=7, ncol=3)
for (i in 6:12) {
  cdc_2017 |>
    filter(Month == i & Experimental.or.control == "Exp.") -> t
  x <- matrix(rep(t$tot.f, 5000), byrow = T, ncol = length(t$tot.f))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,7,T))-mean(x)})
  d_exp[i-5,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$tot.f)
}
plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018),
     xlab="Year", ylab="Population")
lines(2017+(0:6+6)/12, d_con[,2], lwd=2, col=1, lty=2)
arrows(2017+(0:6+6)/12,
       d_con[,1],
       2017+(0:6+6)/12,
       d_con[,3],
       angle=90,
       code=3,
       length=0.05,
       col=1)
lines(2017+(0:6+6)/12, d_exp[,2], lwd=2, col=2, lty=2)
arrows(2017+(0:6+6)/12,
       d_exp[,1],
       2017+(0:6+6)/12,
       d_exp[,3],
       angle=90,
       code=3,
       length=0.05,
       col=2)
grid()
# ----

mali <- MLI
kayes_rural <- single_site(mali, 5)
feed <- seq(0.20,0.24,length.out=5)
out_feed <- list()
for (i in 1:length(feed)) {
  name <- as.character(feed[i])
  kayes_rural_params <- site_parameters(
    interventions = kayes_rural$interventions,
    demography = kayes_rural$demography,
    vectors = kayes_rural$vectors,
    seasonality = kayes_rural$seasonality,
    eir = kayes_rural$eir$eir[1],
    overrides = list(human_population = 5000,
                     mu_atsb = c(feed[i], feed[i], feed[i]))
  )
  kayes_rural_params <- set_atsb(parameters = kayes_rural_params,
                                 timesteps = (17*365+5*30):(18*365), 
                                 coverages = rep(1,366-5*30))
  out_feed[[name]] <- run_simulation(timesteps = kayes_rural_params$timesteps,
                                   parameters = kayes_rural_params)
  print(i)
}

sum_of_squares <- c()
for (i in 1:length(feed)) {
  true_values <- cdc_2017_exp$Mean[4:8]
  y <- out_feed[[i]]$total_M_gambiae*scaler
  x <- out_feed[[i]]$timestep/365 + 2000 + days_shifted/365
  model_values <- approx(x, y, xout= 2017 + c(7:11)/12)
  sum_of_squares[i] <- sum((true_values - model_values$y)^2)
}
plot(feed, sum_of_squares, type="l", lwd=2, frame.plot = F, xlab="Feeding rate")
# minimising sum of squared differences gives excess mortality of 23% as the 
# best fit to the ATSB arm count data

# retrying above using optim function 
out_feed <- list()
sum_squares <- function(feeding_rate) {
  name <- as.character(feeding_rate)
  kayes_rural_params <- site_parameters(
    interventions = kayes_rural$interventions,
    demography = kayes_rural$demography,
    vectors = kayes_rural$vectors,
    seasonality = kayes_rural$seasonality,
    eir = kayes_rural$eir$eir[1],
    overrides = list(human_population = 5000,
                     mu_atsb = c(feeding_rate, feeding_rate, feeding_rate))
    )
  kayes_rural_params <- set_atsb(parameters = kayes_rural_params,
                                 timesteps = (17*365+5*30):(18*365), 
                                 coverages = rep(1,366-5*30))
  out_feed[[name]] <- run_simulation(timesteps = 19*365,
                                     parameters = kayes_rural_params)
  true_values <- cdc_2017_exp$Mean[4:8]
  y <- out_feed[[name]]$total_M_gambiae*scaler
  x <- out_feed[[name]]$timestep/365 + 2000 + days_shifted/365
  model_values <- approx(x, y, xout= 2017 + c(7:11)/12)
  sum_of_squares <- sum((true_values - model_values$y)^2)
  print("1")
  return(sum_of_squares)
}

res1 <- optim(par=0.15, fn = sum_squares, method = "Brent", lower=0.1, upper=0.3, control = list(maxit=5))

plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018),
     xlab="Year", ylab="Population")
lines(out_feed[[4]]$timestep/365+2000+days_shifted/365,
      out_feed[[4]]$total_M_gambiae*scaler,
      col=2,lwd=2)
lines(2017+(cdc_2017_con$Month)/12, cdc_2017_con$Mean, lwd=2, col=1, lty=2)
lines(2017+(cdc_2017_exp$Month)/12, cdc_2017_exp$Mean, lwd=2, col=2, lty=2)

