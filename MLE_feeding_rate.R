################################################################################
## The purpose of this script is to use learn maximum likelihood estimation   ##
## by using using it to fit the feeding rate parameter to the observed count  ##
## data from Mali in 2017. A similar exercise was carried out in Fraser et al.##
## 2020 using a different model and a different method. Hopefully this will   ##
## also give us an idea of how differently they are estimating catch numbers. ##
################################################################################

library(malariasimulation)
library(dplyr)

# load in the count data
cdc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/CDC-Table 1.csv")

# load control malariasim run
malariasim_control <- readRDS("OneDrive_1_2-13-2023/out_mali.RDS")$Kayes_rural_data

# find scaling factor 
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
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
plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018),
     xlab="Year", ylab="Population")
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(out_bounds[[1]]$total_M_gambiae*scaler, rev(out_bounds[[2]]$total_M_gambiae)*scaler),
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