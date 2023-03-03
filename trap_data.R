################################################################################
## This is a test script to examine the trap data for the mali 2016-17 trial
## and find out how to clean it
################################################################################
library(dplyr)
setwd("~/GitHub/atsb_working_code/")
cdc_2016_inside <- read.csv("atsb_working_code/Data Synopsis  catches 2016/CDC 2016 inside-Table 1.csv")
cdc_2016_periphery <- read.csv("atsb_working_code/Data Synopsis  catches 2016/CDC 2016 periphery-Table 1.csv")
malaise_2016_inside <- read.csv("Data Synopsis  catches 2016/Malaise 2016 inside-Table 1.csv")
malaise_2016_periphery <- read.csv("Data Synopsis  catches 2016/Malaise 2016 periphery-Table 1.csv")
psc_2016 <- read.csv("Data Synopsis  catches 2016/PSC 2016-Table 1.csv")

cdc_2016_inside[,-c(19:21)] |>
  filter(!is.na(year)) -> cdc_2016_inside
cdc_2016_periphery[1:204,-c(19:21)] |>
  filter(!is.na(year)) -> cdc_2016_periphery
malaise_2016_inside[,-c(19:21)] |>
  filter(!is.na(year)) -> malaise_2016_inside
malaise_2016_periphery[1:170,-c(19:21)][-c(171:184),] |>
  filter(!is.na(year)) -> malaise_2016_periphery
psc_2016[,-c(19:21)] |>
  filter(!is.na(year)) -> psc_2016

cdc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/CDC-Table 1.csv")
malaise_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/Malaise-Table 1.csv")
psc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/PSC-Table 1.csv")
psc_2017 |>
  filter(!is.na(Month)) -> psc_2017
hlc_apr_jul_2017 <- read.csv("atsb_working_code/HLC data processed/April_May_June_July-Table 1.csv")[,1:16]
hlc_aug_sep_2017 <- read.csv("atsb_working_code/HLC data processed/August_September-Table 1.csv")[,1:16]
hlc_oct_dec_2017 <- read.csv("atsb_working_code/HLC data processed/Oct_Nvber_Dcber-Table 1.csv")[,1:16]
hlc <- rbind(hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)
rm(hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)
# check data is okay 
table(cdc_2016_inside$Month)
table(cdc_2016_periphery$Month)
table(malaise_2016_inside$Month)
table(cdc_2017$Month)

# visualise
plot(cdc_2016_inside$Month, cdc_2016_inside$tot.f, pch=20, frame.plot = F)
points(cdc_2016_periphery$Month, cdc_2016_periphery$tot.f, pch=20, col=2)
points(malaise_2016_inside$Month, malaise_2016_inside$tot.f, pch = 20, col=3)
points(malaise_2016_periphery$Month, malaise_2016_periphery$tot.f, pch=20, col=4)
points(psc_2016$Month, psc_2016$tot.f, pch=20, col=5)

# malaise 2017
plot(malaise_2017$Month, malaise_2017$tot.f, pch=20, frame.plot = F, col=as.factor(cdc_2017$Experimental.or.control),
     ylim=c(0,1000))
malaise_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> malaise_2017_con
lines(malaise_2017_con$Month, malaise_2017_con$Mean, lwd=2, col=1)
arrows(malaise_2017_con$Month, 
       malaise_2017_con$Mean-malaise_2017_con$CI,
       malaise_2017_con$Month,
       malaise_2017_con$Mean+malaise_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=1)
malaise_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> malaise_2017_exp
lines(malaise_2017_exp$Month, malaise_2017_exp$Mean, lwd=2, col=2)
arrows(malaise_2017_exp$Month, 
       malaise_2017_exp$Mean-malaise_2017_exp$CI,
       malaise_2017_exp$Month,
       malaise_2017_exp$Mean+malaise_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=2)

# psc 2017
plot(psc_2017$Month, psc_2017$tot.f, pch=20, frame.plot = F, col=as.factor(cdc_2017$Experimental.or.control))
psc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> psc_2017_con
lines(psc_2017_con$Month, psc_2017_con$Mean, lwd=2, col=1)
arrows(psc_2017_con$Month, 
       psc_2017_con$Mean-psc_2017_con$CI,
       psc_2017_con$Month,
       psc_2017_con$Mean+psc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=1)
psc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> psc_2017_exp
lines(psc_2017_exp$Month, psc_2017_exp$Mean, lwd=2, col=2)
arrows(psc_2017_exp$Month, 
       psc_2017_exp$Mean-psc_2017_exp$CI,
       psc_2017_exp$Month,
       psc_2017_exp$Mean+psc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=2)

# cdc 2017
plot(cdc_2017$Month, cdc_2017$tot.f, pch=20, frame.plot = F, col=as.factor(cdc_2017$Experimental.or.control))
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> cdc_2017_con
lines(cdc_2017_con$Month, cdc_2017_con$Mean, lwd=2, col=1)
arrows(cdc_2017_con$Month, 
       cdc_2017_con$Mean-cdc_2017_con$CI,
       cdc_2017_con$Month,
       cdc_2017_con$Mean+cdc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=1)
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> cdc_2017_exp
lines(cdc_2017_exp$Month, cdc_2017_exp$Mean, lwd=2, col=2)
arrows(cdc_2017_exp$Month, 
       cdc_2017_exp$Mean-cdc_2017_exp$CI,
       cdc_2017_exp$Month,
       cdc_2017_exp$Mean+cdc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=2)

cdc_2017 |>
  filter(Month > 5) |>
  mutate(Date=paste0(Date, "_", Month, "_2016")) |>
  mutate(Date=parse_date_time(Date, "dmy")) -> cdc_2017
plot(rep(2017,nrow(cdc_2017)), cdc_2017$tot.f, pch=20, frame.plot = F, col=as.factor(cdc_2017$Experimental.or.control))
cdc_2017 |>
  group_by(Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> cdc_2017_con
arrows(2017, 
       cdc_2017_con$Mean-cdc_2017_con$CI,
       2017,
       cdc_2017_con$Mean+cdc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=1)
cdc_2017 |>
  group_by(Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> cdc_2017_exp
arrows(2017, 
       cdc_2017_exp$Mean-cdc_2017_exp$CI,
       2017,
       cdc_2017_exp$Mean+cdc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=2)
points(malaise_2017$Month, malaise_2017$tot.f, pch=20, col=as.factor(malaise_2017$Experimental.or.control))
points(psc_2017$Month, psc_2017$tot.f, pch=20, col=as.factor(psc_2017$Experimental.or.control))

# hlc 
plot(hlc_2017$Month, hlc_2017$Count, pch=20, frame.plot = F, col=as.factor(hlc_2017$Treatment))
hlc_2017 |>
  group_by(Month, Treatment) |>
  summarise(Mean=mean(Count),CI=1.96*sd(Count)/sqrt(n())) |>
  filter(Treatment==0) -> hlc_2017_con
lines(hlc_2017_con$Month, hlc_2017_con$Mean, lwd=2, col=1)
arrows(hlc_2017_con$Month, 
       hlc_2017_con$Mean-hlc_2017_con$CI,
       hlc_2017_con$Month,
       hlc_2017_con$Mean+hlc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=1)
hlc_2017 |>
  group_by(Month, Treatment) |>
  summarise(Mean=mean(Count),CI=1.96*sd(Count)/sqrt(n())) |>
  filter(Treatment==1) -> hlc_2017_exp
lines(hlc_2017_exp$Month, hlc_2017_exp$Mean, lwd=2, col=2)
arrows(hlc_2017_exp$Month, 
       hlc_2017_exp$Mean-hlc_2017_exp$CI,
       hlc_2017_exp$Month,
       hlc_2017_exp$Mean+hlc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=2)
title("Human landing catch")

# glms
cdc_2017$Vilage <- as.factor(cdc_2017$Vilage)
cdc_2017$Date <- as.factor(cdc_2017$Date)
cdc_2017$Month <- as.factor(cdc_2017$Month)
fit <- glmer.nb(tot.f ~ as.factor(Experimental.or.control) + (1|Date), data = cdc_2017) 
summary(fit)
exp(2.6698+3.44-1.205)
rho1_2 <- vcov(fit)[1,2]/(sqrt(vcov(fit)[1,1])*sqrt(vcov(fit)[2,2]))
sigma1_2 <-  sqrt(vcov(fit)[1,1] + vcov(fit)[2,2] + 2 * rho1_2 *(sqrt(vcov(fit)[1,1]) *(sqrt(vcov(fit)[2,2]))))


glmer.nb(tot.f ~ as.factor(Experimental.or.control) + (1|Vilage) + (1|Month), data = cdc_2017) |>
  summary()
exp(5.2250)

fit <- glm.nb(tot.f ~ as.factor(Experimental.or.control), data = malaise_2017) 
summary(fit)
exp(4.8845)
