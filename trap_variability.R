setwd("Documents/GitHub/")

library(lme4)
library(dplyr)

# load in trap data here
cdc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/CDC-Table 1.csv")

cdc_2017_bytrap <- matrix(0, nrow=1260, ncol=6)

for (i in 1:nrow(cdc_2017)) {
  for (j in 1:10) {
    cdc_2017_bytrap[(i-1)*10+j,6] <- cdc_2017[i,j+7]
    cdc_2017_bytrap[(i-1)*10+j,5] <- j
    cdc_2017_bytrap[(i-1)*10+j,4] <- cdc_2017[i,4]
    cdc_2017_bytrap[(i-1)*10+j,3] <- cdc_2017[i,3]
    cdc_2017_bytrap[(i-1)*10+j,2] <- cdc_2017[i,2]
    cdc_2017_bytrap[(i-1)*10+j,1] <- cdc_2017[i,1]
  }
}
colnames(cdc_2017_bytrap) <- c("Month", "Date", "Vilage", "Experimental.or.control",
                               "Trap_number", "Count")
cdc_2017_bytrap <- data.frame(cdc_2017_bytrap)
cdc_2017_bytrap$Observation <- as.factor(1:1260)
cdc_2017_bytrap$Month <- as.numeric(cdc_2017_bytrap$Month)
cdc_2017_bytrap$Count <- as.numeric(cdc_2017_bytrap$Count)

# visualise
plot(cdc_2017_bytrap$Month, cdc_2017_bytrap$Count, col=as.factor(cdc_2017_bytrap$Experimental.or.control),
     frame.plot = F, pch=20, cex=0.7, xlab = "Month", ylab = "Count")
cdc_2017_bytrap |>
  filter(Experimental.or.control=="Con.") |>
  group_by(Month) |>
  summarise(Mean=mean(Count), CI=1.96*sd(Count)/sqrt(n())) -> cdc_2017_bytrap_con
cdc_2017_bytrap |>
  filter(Experimental.or.control=="Exp.") |>
  group_by(Month) |>
  summarise(Mean=mean(Count), CI=1.96*sd(Count)/sqrt(n())) -> cdc_2017_bytrap_exp

lines(cdc_2017_bytrap_con$Month, cdc_2017_bytrap_con$Mean, lwd=2, col=1)
lines(cdc_2017_bytrap_exp$Month, cdc_2017_bytrap_exp$Mean, lwd=2, col=2)
arrows(cdc_2017_bytrap_con$Month, 
       cdc_2017_bytrap_con$Mean-cdc_2017_bytrap_con$CI,
       cdc_2017_bytrap_con$Month,
       cdc_2017_bytrap_con$Mean+cdc_2017_bytrap_con$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=1,
       lwd=2)
arrows(cdc_2017_bytrap_exp$Month, 
       cdc_2017_bytrap_exp$Mean-cdc_2017_bytrap_exp$CI,
       cdc_2017_bytrap_exp$Month,
       cdc_2017_bytrap_exp$Mean+cdc_2017_bytrap_exp$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col=2,
       lwd=2)

#### 
d_con <- matrix(0, nrow=7, ncol=3)
for (i in 6:12) {
  cdc_2017_bytrap |>
    filter(Month == i & Experimental.or.control == "Con.") -> t
  x <- matrix(rep(t$Count, 5000), byrow = T, ncol = length(t$Count))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,70,T))-mean(x)})
  d_con[i-5,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$Count)
}

d_exp <- matrix(0, nrow=7, ncol=3)
for (i in 6:12) {
  cdc_2017_bytrap |>
    filter(Month == i & Experimental.or.control == "Exp.") -> t
  x <- matrix(rep(t$Count, 5000), byrow = T, ncol = length(t$Count))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,70,T))-mean(x)})
  d_exp[i-5,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$Count)
}
plot(cdc_2017_bytrap$Month, cdc_2017_bytrap$Count, col=as.factor(cdc_2017_bytrap$Experimental.or.control),
     frame.plot = F, pch=20, cex=0.6, xlab = "Month", ylab = "Count")
lines(6:12, d_con[,2], lwd=2, col=1, lty=1)
arrows(6:12,
       d_con[,1],
       6:12,
       d_con[,3],
       angle=90,
       code=3,
       length=0.1,
       col=1,
       lwd=2)
lines(6:12, d_exp[,2], lwd=2, col=2, lty=1)
arrows(6:12,
       d_exp[,1],
       6:12,
       d_exp[,3],
       angle=90,
       code=3,
       length=0.1,
       col=2,
       lwd=2)
legend(x="topright", legend=c("Control", "ATSB"), col=c(1,2), lty=1, bty="n", lwd=2)

cdc_2017_bytrap |>
  filter(Month > 6 & Month < 12) -> cdc_2017_bytrap

cdc_2017_bytrap$Month <- as.factor(cdc_2017_bytrap$Month)
cdc_2017_bytrap$Vilage <- as.factor(cdc_2017_bytrap$Vilage)
cdc_2017_bytrap$Trap_number <- as.factor(cdc_2017_bytrap$Trap_number)
cdc_2017_bytrap$Count <- as.numeric(cdc_2017_bytrap$Count)

fit <- glmer.nb(Count ~ as.factor(Experimental.or.control) + Month 
           + (1|Vilage) + (1|Trap_number) + (1|Observation), data=cdc_2017_bytrap) |> summary()
rho1_2 <- vcov(fit)[1,2]/(sqrt(vcov(fit)[1,1])*sqrt(vcov(fit)[2,2]))
sigma1_2 <-  sqrt(vcov(fit)[1,1] + vcov(fit)[2,2] + 2 * rho1_2 *(sqrt(vcov(fit)[1,1]) *(sqrt(vcov(fit)[2,2]))))
summary(fit)

cdc_2017_bytrap |> 
  filter(Experimental.or.control=="Exp.") |>
  summarise(Mean=mean(Count))





 