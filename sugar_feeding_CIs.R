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

# read in the sugar feeding data
sugar_feeding <- read.csv("DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
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
  filter(month >= 6) -> sugar_feeding 

sugar_feeding |>
  group_by(month) |>
  summarise(feeding_rate=females.ASB.positive/TOTAL.Sample.females.Day.2) -> t

x <- matrix(rep(t$feeding_rate,5000), byrow=T, ncol=length(t$feeding_rate)) # these two lines took me 2 hours 
a <- apply(x, MARGIN=1, FUN = function(x){mean(sample(x,49,T))-mean(x)})    # for some reason, i was trying to learn
quantile(a, c(0.025,0.975))                                                 # to use sapply but then realised I could use apply


data <- t$females.ASB.positive/t$TOTAL.Sample.females.Day.2
n <- 5000
d <- numeric(n)
for (i in 1:n) {
 s <- sample(data, 3, replace = TRUE)
 d[i] <- mean(s)-mean(data)
}
hist(d, breaks = 30)
quantile(d, probs=c(0.025,0.975))

x <- matrix(rep(data,5000), byrow=T, ncol=length(data))
a <- apply(x, MARGIN=1, FUN = function(x){mean(sample(x,7,T))-mean(x)})
quantile(a, c(0.025,0.975))
