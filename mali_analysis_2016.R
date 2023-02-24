library(tidyverse)
library(lubridate)
require(pscl)
require(MASS)
require(boot)
library(lme4)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(fitdistrplus)
InvLogit <- function(X){
  exp(X)/(1+exp(X))
}

setwd("Documents/GitHub/")
hlc_may_jul <- read.csv("HLC data 2016 (unchanged) processed/HLC May_Jul2016-Table 1.csv") 
hlc_aug_sep <- read.csv("HLC data 2016 (unchanged) processed/HLC Aug_Sept 2016-Table 1.csv")
hlc_oct_dec <- read.csv("HLC data 2016 (unchanged) processed/HLC Oct_Dber-Table 1.csv")
hlc_aug_sep <- hlc_aug_sep[,-17]
# the data for 2017 has no date so havent included it yet
#hlc_apr_jul_2017 <- read.csv("HLC data processed/April_May_June_July-Table 1.csv")[,1:16]
# hlc_aug_sep_2017 <- read.csv("HLC data processed/August_September-Table 1.csv")[,1:16]
# hlc_oct_dec_2017 <- read.csv("HLC data processed/Oct_Nvber_Dcber-Table 1.csv")[,1:16]

hlc <- rbind(hlc_may_jul, hlc_aug_sep, hlc_oct_dec)
rm(hlc_may_jul, hlc_aug_sep, hlc_oct_dec)

hlc <- hlc[-which(hlc$NightFrames=="N/A"),]
hlc$Date[1:63] <- c("1-May-2016")
hlc$Date[64:181] <- c("1-June-2016")

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
  summarise(Count=n(), Outdoor=sum(Position=="Outside"), Indoor=sum(Position=="Inside")) -> hlc_bysitedate

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
  group_by(Date) |>
  summarise(Count=n(), Outdoor=sum(Position=="Outside"), Indoor=sum(Position=="Inside")) -> hlc_bydate

hlc_bysitedate |>
  mutate(Date=parse_date_time(Date, "dmy")) |>
  arrange(Date) |>
  mutate(Month_Yr = format_ISO8601(Date, precision = "ym")) -> hlc_bysitedate
hlc_bysitedate$Site <- as.factor(hlc_bysitedate$Site)
hlc_bysitedate$Month_Yr <- as.factor(hlc_bysitedate$Month_Yr)
hlc_bysitedate$Date <- as.factor(hlc_bysitedate$Date)

hlc_bydate |>
  mutate(Date=parse_date_time(Date, "dmy")) |>
  arrange(Date) |>
  mutate(Month_Yr = format_ISO8601(Date, precision = "ym")) -> hlc_bydate
hlc_bydate$Month_Yr <- as.factor(hlc_bydate$Month_Yr)
hlc_bydate$Date <- as.factor(hlc_bydate$Date)

# visualising data
plot.default(hlc_bydate$Month_Yr, hlc_bydate$Count, pch=20, frame.plot=F)
hlc_bydate |> 
  group_by(Month_Yr) |>
  summarise(mean=mean(Count), sd=sd(Count)) -> out
lines(1:8, out$mean, col="red", lwd=2)
arrows(1:8, 
       out$mean-out$sd,
       1:8,
       out$mean+out$sd,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")



# Mixed and random effects regressions
# NULL model - no fixed effect and random effect on date
glmer.nb(Count ~ (1|Date), data = hlc_bydate) |> 
  summary()
glmer.nb(Count ~ (1|Date), data = hlc_bysitedate) |> 
  summary()

glmer(
  cbind(Indoor, Count - Indoor) ~ (1 | Date),
  family = binomial, data = hlc_bydate) |>
  summary()
glmer(
  cbind(Indoor, Count - Indoor) ~ (1 | Date),
  family = binomial, data = hlc_bysitedate) |>
  summary()

# NULL model - no fixed effect and random effects on date and site
glmer.nb(Count ~ (1|Date) + (1|Site), data = hlc_bysitedate) |>
  summary()
glmer(
  cbind(Indoor, Count - Indoor) ~ (1 | Date) + (1|Site),
  family = binomial, data = hlc_bysitedate) |>
  summary()

# fixed effect on month random effect on date
glmer.nb(Count ~ Month_Yr + (1|Date), data = hlc_bydate) |>
  summary()
mixed_model(fixed = Count ~ Month_Yr, random = ~ 1|Date, data=hlc_bysitedate,
            family = negative.binomial()) |> summary()

glmer(
  cbind(Indoor, Count - Indoor) ~ Month_Yr + (1 | Date),
  family = binomial, data = hlc_bydate) |>
  summary()
glmer(
  cbind(Indoor, Count - Indoor) ~ Month_Yr + (1 | Date),
  family = binomial, data = hlc_bysitedate) |>
  summary()

# fixed effect on month random effects on date and site
glmer.nb(Count ~ Month_Yr + (1|Date) + (1|Site), data = hlc_bysitedate) |>
  summary()


sugar_feeding <- read.csv("DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
sugar_feeding <- sugar_feeding[-(64:65),]
table(sugar_feeding$Village)
mean(sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2)
sd(sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2)

par(las=1)
plot.default(sugar_feeding$month, 
             sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2,
             cex=0.9, pch=20, frame.plot = F, xlab = "Month", ylab = "% bait fed",
             ylim = c(0,1))
qt(c(0.025,0.975), df = 6)
sugar_feeding |>
  group_by(month) |>
  summarise(mean = mean(females.ASB.positive/TOTAL.Sample.females.Day.2),
             CI = 2.44*sd(females.ASB.positive/TOTAL.Sample.females.Day.2)/sqrt(n())) -> sugar_feeding_grouped

lines(sugar_feeding_grouped$month, sugar_feeding_grouped$mean, col="red", lwd=2)
arrows(sugar_feeding_grouped$month, 
       sugar_feeding_grouped$mean-sugar_feeding_grouped$CI,
       sugar_feeding_grouped$month,
       sugar_feeding_grouped$mean+sugar_feeding_grouped$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")

plot.default(sugar_feeding$Village, 
             sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2,
             cex=0.9, pch=20, frame.plot = F, xlab = "Village", ylab = "% bait fed",
             ylim = c(0,1))
qt(c(0.025,0.975), df = 8)
sugar_feeding |>
  group_by(Village) |>
  summarise(mean = mean(females.ASB.positive/TOTAL.Sample.females.Day.2),
            CI = 2.31*sd(females.ASB.positive/TOTAL.Sample.females.Day.2)/sqrt(n())) -> sugar_feeding_grouped
lines(sugar_feeding_grouped$Village, sugar_feeding_grouped$mean, col="red", lwd=2)
arrows(1:7, 
       sugar_feeding_grouped$mean-sugar_feeding_grouped$CI,
       1:7,
       sugar_feeding_grouped$mean+sugar_feeding_grouped$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")

plot(rep(1,63), 
     sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2,
     cex=0.9, pch=20, frame.plot = F, xlab = "Month", ylab = "% bait fed",
     ylim = c(0,1))
qt(c(0.025,0.975), df = 62)
sugar_feeding |>
  summarise(mean = mean(females.ASB.positive/TOTAL.Sample.females.Day.2),
            CI = 2.00*sd(females.ASB.positive/TOTAL.Sample.females.Day.2)/sqrt(n())) -> sugar_feeding_grouped
lines(1, sugar_feeding_grouped$mean, col="red", lwd=2)
arrows(1, 
       sugar_feeding_grouped$mean-sugar_feeding_grouped$CI,
       1,
       sugar_feeding_grouped$mean+sugar_feeding_grouped$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")
sugar_feeding |>
  filter(month>5) -> sugar_feeding
sugar_feeding$month <- as.factor(sugar_feeding$month)
sugar_feeding$Village <- as.factor(sugar_feeding$Village)

fit <- glmer(
  cbind(females.ASB.positive, TOTAL.Sample.females.Day.2 - females.ASB.positive) ~ month + (1 | Village),
  family = binomial, data = sugar_feeding) 
summary(fit)
glm(
  cbind(females.ASB.positive, TOTAL.Sample.females.Day.2 - females.ASB.positive) ~ month,
      family = binomial, data = sugar_feeding) |> 
      summary()
InvLogit(-0.10996+0.25*c(1.96, 0, -1.96))
InvLogit(-0.1694-0.2221+0.3349*c(1.96, 0, -1.96))
InvLogit(-0.1694-0.7953+0.2432*c(1.96, 0, -1.96))
