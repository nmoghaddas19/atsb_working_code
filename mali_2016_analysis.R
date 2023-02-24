library(tidyverse)
library(lubridate)
setwd("Documents/GitHub/")

hlc_may_jul <- read.csv("HLC data 2016 (unchanged) processed/HLC May_Jul2016-Table 1.csv") 
hlc_aug_sep <- read.csv("HLC data 2016 (unchanged) processed/HLC Aug_Sept 2016-Table 1.csv")
hlc_oct_dec <- read.csv("HLC data 2016 (unchanged) processed/HLC Oct_Dber-Table 1.csv")
hlc_aug_sep <- hlc_aug_sep[,-17]
hlc_apr_jul_2017 <- read.csv("HLC data processed/April_May_June_July-Table 1.csv")[,1:16]
hlc_aug_sep_2017 <- read.csv("HLC data processed/August_September-Table 1.csv")[,1:16]
hlc_oct_dec_2017 <- read.csv("HLC data processed/Oct_Nvber_Dcber-Table 1.csv")[,1:16]

colnames(hlc_apr_jul_2017) <- colnames(hlc_may_jul)
colnames(hlc_aug_sep_2017) <- colnames(hlc_may_jul)
colnames(hlc_oct_dec_2017) <- colnames(hlc_may_jul)
hlc <- rbind(hlc_may_jul, hlc_aug_sep, hlc_oct_dec)

rm(hlc_may_jul, hlc_aug_sep, hlc_oct_dec)

unique(hlc$Study.Sites)
unique(hlc$Month)
unique(hlc$Date)
unique(hlc$Species)
unique(hlc$NightFrames)
unique(hlc$Parity)
unique(hlc$fed.unfed)
sum(hlc$fed == "N/A")
sum(hlc$NightFrames=="N/A")

hlc <- hlc[-which(hlc$NightFrames=="N/A"),]
hlc$Date[1:63] <- c("1-May-2016")
hlc$Date[64:181] <- c("1-June-2016")

table(hlc$Position)

hlc |>
  mutate(time_frame = ifelse(NightFrames=="18h-20h", 1,
                             ifelse(NightFrames=="20-22h"|NightFrames=="20h-22h", 2,
                                    ifelse(NightFrames=="22h-00h"|NightFrames=="22H-00h", 3,
                                           ifelse(NightFrames=="00h-02h", 4,
                                                  ifelse(NightFrames=="02h-04h"|NightFrames=="02h-4h", 5,
                                                         ifelse(NightFrames=="04h-06h"|NightFrames=="04H-06h", 6,
                                                                NightFrames))))))) |> 
  mutate(Site = ifelse(Study.Sites=="Kignélé"|Study.Sites=="Kignele ", "Kignele",
                       ifelse(Study.Sites=="Niaganabougou"|Study.Sites=="Nianganabougou", "Nianguanabougou",
                              ifelse(Study.Sites=="Sirakélé","Sirakele",
                                     ifelse(Study.Sites=="Krékrélo","Krekrelo",
                                            ifelse(Study.Sites=="Trekourou", "Trekrou",
                                                   Study.Sites)))))) |>
  mutate(Parity = ifelse(Parity=="p"|Parity=="PP", "P",
                         ifelse(Parity=="NO"|Parity=="MP", "NP",
                                Parity))) |>
  group_by(Site, Date) |>
  summarise(Count=n(), Outdoor=sum(Position=="Outside"), Indoor=sum(Position=="Inside")) -> hlc

table(hlc$time_frame)
table(hlc$Study.Sites)
table(hlc$Month)
table(hlc$fed.unfed)

hlc %>%
  group_by(Site, Date, time_frame) %>%
  summarise(Hourly_Count = sum(Count)) %>%
  ungroup() %>%
  complete(Site, Date, time_frame, fill = list(Hourly_Count = 0)) -> hlc

stupid_date <- parse_date_time(Date, "dmy")
hlc |>
  mutate(Date=parse_date_time(Date, "dmy")) |>
  arrange(Date) |>
  mutate(Month_Yr = format_ISO8601(Date, precision = "ym")) -> hlc_bysitedate

hlc$Site <- as.factor(hlc$Site)
hlc$Month_Yr <- as.factor(hlc$Month_Yr)
fig.pois <- fitdist(hlc$Count, "pois") # Poisson distr
plot(fig.pois)
fig.negbin <- fitdist(hlc$Count, "nbinom") # negative binomial
plot(fig.negbin)

plot.default(hlc_bysitedate$Month_Yr, hlc_bysitedate$Count, pch=20, cex=0.9, frame.plot = F,
             xaxt = "n", xlab = "Month", ylab = "Daily Count")
axis(1, at=1:8, labels=c("May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"))
hlc_bysitedate |> 
  group_by(Month_Yr) |>
  summarise(mean=mean(Count), sd=sd(Count)) -> outt
lines(1:8, outt$mean, col="red", lwd=2)
arrows(1:8, 
       outt$mean-outt$sd,
       1:8,
       outt$mean+outt$sd,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")

fit_nb_bysitedate <- glmer.nb(Count ~ Month_Yr + (1|Site), data = hlc_bysitedate)
summary(fit_nb_bysitedate)
exp(coef(summary(fit_nb))["(Intercept)", "Estimate"] + coef(summary(fit_nb))["Month_Yr2016-08", "Estimate"])
rho1_3 <- vcov(fit_nb)[1,3]/(sqrt(vcov(fit_nb)[1,1])*sqrt(vcov(fit_nb)[3,3]))
sigma <- sqrt(vcov(fit_nb)[1,1] + vcov(fit_nb)[3,3] + 2 * rho1_3 *(sqrt(vcov(fit_nb)[1,1]) *(sqrt(vcov(fit_nb)[3,3]))))
require(ggplot2)
install.packages("pscl")
require(pscl)
require(MASS)
require(boot)

hlc$time_frame <- as.numeric(hlc$time_frame)

table(hlc$Month_Yr)/6
table(completed.df_1$Month_Yr)/6

graph(hlc)
title("Hourly Counts")
months_2016 <- unique(hlc$Month_Yr)

plot.new()
par(mfrow=c(3,3))
for (i in 1:length(months_2016)) {
  hlc |>
    filter(grepl(months_2016[i], Date)) |>
    graph()
  title(months_[i])
}

hlc |>
  group_by(Date) |>
  summarise(Daily_Count = sum(Hourly_Count)) |> 
  mutate(Month_Yr = as.factor(format_ISO8601(Date, precision = "ym"))) -> Daily_Counts_2016

plot.new()
par(mfrow=c(1,1))
hist(Daily_Counts_2016$Daily_Count, breaks=15, xlab = "Daily Count")

plot.default(Daily_Counts_2016$Month_Yr, Daily_Counts_2016$Daily_Count, pch=20, cex=0.9, frame.plot = F,
             xaxt = "n", xlab = "Month", ylab = "Daily Count")
axis(1, at=1:8, labels=months_2016)
Daily_Counts_2016 |> 
  group_by(Month_Yr) |>
  summarise(mean=mean(Daily_Count), sd=sd(Daily_Count)) -> outt
lines(1:8, outt$mean, col="red", lwd=2)
arrows(1:8, 
       outt$mean-outt$sd,
       1:8,
       outt$mean+outt$sd,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")
title("Daily Counts by Month")

library(lme4)
library(MASS)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(fitdistrplus)

fig.pois <- fitdist(hlc_bydate$Count, "pois") # Poisson distr
plot(fig.pois)
fig.negbin <- fitdist(hlc_bydate$Count, "nbinom") # negative binomial
plot(fig.negbin)
#Here we see the Negative binomial model has a lower Akaike's Information Criterion,
# Hence a better fit
gofstat(list(fig.pois,fig.negbin), fitnames = c("Poisson", "Negative Binomial"))

library(GLMMmisc)
InvLogit <- function(X){
  exp(X)/(1+exp(X))
}
install.packages("GLMMadaptive")
library(GLMMadaptive)

hlc_bydate$Date <- as.factor(hlc_bydate$Date)
hlc_bydate$Month_Yr <- as.factor(hlc_bydate$Month_Yr)
hlc$Month_Yr <- as.factor(hlc$Month_Yr)
hlc$Date <- as.factor(hlc$Date)
fit_nb_bydate <- glmer.nb(Count ~ Month_Yr + (1|Day), data = hlc_bydate)
summary(fit_nb_bydate)
exp(coef(summary(fit_nb_bydate))["(Intercept)", "Estimate"] + coef(summary(fit_nb_bydate))["Month_Yr2016-08", "Estimate"])
rho1_3 <- vcov(fit_nb_bydate)[1,3]/(sqrt(vcov(fit_nb_bydate)[1,1])*sqrt(vcov(fit_nb_bydate)[3,3]))
sigma <- sqrt(vcov(fit_nb_bydate)[1,1] + vcov(fit_nb_bydate)[3,3] + 2 * rho1_3 *(sqrt(vcov(fit_nb_bydate)[1,1]) *(sqrt(vcov(fit_nb_bydate)[3,3]))))
exp(0.8079)

?mixed_model
?glmer.nb
mixed_model(fixed = Count ~ Month_Yr, random =  ~ 1|Date, data=hlc_bydate, 
             family = negative.binomial()) |> summary()
glmer.nb(Count ~ Month_Yr + (1|Date), data = hlc_bydate) |> 
  summary()
mixed_model(fixed = Count ~ Site, random = ~ 1|Date, data=hlc_bysitedate,
            family = negative.binomial()) |> summary()
glmer.nb(Count ~ Site + (1|Date), data = hlc_bysitedate) |> 
  summary()

glmer(
  cbind(Indoor, Count - Indoor) ~ 
   (1 | Date),
  family = binomial, data = hlc_bydate) |>
  summary()

InvLogit <- function(X){
  exp(X)/(1+exp(X))
}

InvLogit(1.584)

Daily_Counts$Date <- as.factor(Daily_Counts$Date)
fit_nb <- glmer.nb(Daily_Count ~ Month_Yr + (1|Date), data = Daily_Counts)
summary(fit_nb)

