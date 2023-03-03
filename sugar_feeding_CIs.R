################################################################################
## The purpose of this script is to generate confidence intervals for the     ##
## ASB sugar feeding data. I will attempt to do this by fitting a mixed       ##      
## effects GLM and also using a bootstrap method. The aim is to use the       ##
## confidence intervals to generate predictions in mosquito count numbers     ##
## from malariasimulation to which I can compare the measured count data from ##
## the 2016-2017 Mali entomological trial.                                    ##
################################################################################

# packages
library(lubridate)
require(pscl)
require(MASS)
require(boot)
library(lme4)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(fitdistrplus)
library(foresite)
library(site)

# read in the sugar feeding data
sugar_feeding <- read.csv("atsb_working_code/DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
sugar_feeding <- sugar_feeding[-(64:65),] # final two rows are NAs so removing them

# visualise the data
par(las=1)
table(sugar_feeding$month)
plot.default(sugar_feeding$month, 
             sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2,
             cex=0.9, pch=20, frame.plot = F, xlab = "Month", ylab = "% bait fed",
             ylim = c(0,1))
qt(c(0.025,0.975), df = 6) # t distribution multiplier for 95% CIs with 6 degrees of freedom
sugar_feeding |>
  group_by(month) |>
  summarise(mean = mean(females.ASB.positive/TOTAL.Sample.females.Day.2),
            CI = 2.44*sd(females.ASB.positive/TOTAL.Sample.females.Day.2)/sqrt(n())) -> sugar_feeding_grouped
polygon(c(4,12,12,4), c(0.15,0.15,0.45,0.45),
        col = adjustcolor("grey", alpha.f = 0.5), border = F)
lines(sugar_feeding_grouped$month, sugar_feeding_grouped$mean, col="red", lwd=2)
arrows(sugar_feeding_grouped$month, 
       sugar_feeding_grouped$mean-sugar_feeding_grouped$CI,
       sugar_feeding_grouped$month,
       sugar_feeding_grouped$mean+sugar_feeding_grouped$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")

# lets filter the data from month 6-12 since thats when the ATSBs were active
sugar_feeding |>
  filter(month > 5) |>
  mutate(observation=as.factor(1:49)) -> sugar_feeding 
sugar_feeding$month <- as.factor(sugar_feeding$month)
sugar_feeding$Village <- as.factor(sugar_feeding$Village)

# alright lets try fitting some glms to estimate confidence intervals
glmer(
  cbind(females.ASB.positive, TOTAL.Sample.females.Day.2 - females.ASB.positive) ~ (1|Village) + (1|month) + (1|observation),
  family = binomial, data = sugar_feeding) |>
  summary()
InvLogit(-0.55492+0.02)
InvLogit(-0.1694-0.2221+0.3349*c(1.96, 0, -1.96))

# alright lets use bootstrap to generate some confidence intervals
# grouping months together
sugar_feeding |>
  group_by(month) |>
  summarise(feeding_rate=females.ASB.positive/TOTAL.Sample.females.Day.2) -> t

x <- matrix(rep(t$feeding_rate,5000), byrow=T, ncol=length(t$feeding_rate))      # these two lines took me 2 hours 
a <- apply(x, MARGIN=1, FUN = function(x){mean(sample(x,nrow(t),T))-mean(x)})    # for some reason, I was trying to learn
quantile(a, c(0.025,0.975))                                                      # to use sapply but then realised I could use apply
bounds <- c(mean(t$feeding_rate), mean(t$feeding_rate)+quantile(a, c(0.025,0.975)))
bounds <- c(0.078,0.38)
bounds <- c(0.27, 0.395)

# disaggregating by month
sugar_feeding_CI <- matrix(0,nrow=7,ncol=3)
for (i in 6:12) {
  sugar_feeding |>
    filter(month==i) |>
    group_by(month) |>
    summarise(feeding_rate=females.ASB.positive/TOTAL.Sample.females.Day.2) -> this_sugar_feeding
  x <- matrix(rep(this_sugar_feeding$feeding_rate,5000), byrow=T, 
              ncol=length(this_sugar_feeding$feeding_rate)) 
  a <- apply(x, MARGIN=1, FUN = function(x){mean(sample(x,nrow(this_sugar_feeding),T))-mean(x)})
  sugar_feeding_CI[i-5,] <- c(mean(this_sugar_feeding$feeding_rate), 
                              mean(this_sugar_feeding$feeding_rate)+quantile(a, c(0.025,0.975)))
}

matrix(c(0.6047629, 0.3342333, 0.1414194,
         0.4703416, 0.3452918, 0.2385171,
         0.3913919, 0.2734477, 0.1805043,
         0.5136576, 0.3460608, 0.2095825,
         0.4378784, 0.2979038, 0.1877307,
         0.5081354, 0.3717305, 0.2531004,
         0.6122226, 0.3614214, 0.1686725), byrow=T, ncol=3) |>
  apply(MARGIN=2, FUN=mean) -> monthly
bounds <- monthly[c(1,3)]

# okay what range of counts does malariasim predict for those confidence intervals?
mali <- MLI
kayes_rural <- single_site(mali, 5)
out_bounds <- list()
for (i in 1:length(bounds)) {
  name <- as.character(bounds[i])
  kayes_rural_params <- site_parameters(
    interventions = kayes_rural$interventions,
    demography = kayes_rural$demography,
    vectors = kayes_rural$vectors,
    seasonality = kayes_rural$seasonality,
    eir = kayes_rural$eir$eir[1],
    overrides = list(human_population = 5000,
                     mu_atsb = c(bounds[i], bounds[i], bounds[i]))
  )
  kayes_rural_params <- set_atsb(parameters = kayes_rural_params,
                                 timesteps = (17*365+5*30):(18*365), 
                                 coverages = rep(1,366-5*30))
  out_bounds[[name]] <- run_simulation(timesteps = kayes_rural_params$timesteps,
                                      parameters = kayes_rural_params)
  print(i)
}

par(las=1)
plot(o[[1]]$timestep/365+2000,
     o[[1]]$total_M_gambiae/180.0287, type="l", lwd=2, frame.plot = F, xlim=c(2016,2018),
     xlab="Year", ylab="Count", ylim=c(0,800))
polygon(c(out_bounds[[1]]$timestep/365+2000, rev(out_bounds[[2]]$timestep/365+2000)),
        c(out_bounds[[1]]$total_M_gambiae/180.0287, rev(out_bounds[[2]]$total_M_gambiae)/180.0287),
        col = adjustcolor("dodgerblue", alpha.f = 0.4), border = FALSE)

# points(2017+cdc_2017$Month/12, cdc_2017$tot.f, pch=20, col=as.factor(cdc_2017$Experimental.or.control))
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> cdc_2017_con
lines(2017+(cdc_2017_con$Month-1)/12, cdc_2017_con$Mean, lwd=2, col=1, lty=2)
arrows(2017+(cdc_2017_con$Month-1)/12,
       cdc_2017_con$Mean-cdc_2017_con$CI,
       2017+(cdc_2017_con$Month-1)/12,
       cdc_2017_con$Mean+cdc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.05,
       col=1)
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> cdc_2017_exp
lines(2017+(cdc_2017_exp$Month-1)/12, cdc_2017_exp$Mean, lwd=2, col=2, lty=2)
arrows(2017+(cdc_2017_exp$Month-1)/12, 
       cdc_2017_exp$Mean-cdc_2017_exp$CI,
       2017+(cdc_2017_exp$Month-1)/12,
       cdc_2017_exp$Mean+cdc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=2)
legend(x="topright", legend=c("Control", "ATSB"), col=c(1,2), lwd=2, lty=2, bty="n")


