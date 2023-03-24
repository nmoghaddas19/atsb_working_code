# Looking at anthrone feeding and if it correlates to ASB feeding

library(dplyr)
library(boot)

# read in the sugar feeding data
sugar_feeding <- read.csv("atsb_working_code/DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
sugar_feeding <- sugar_feeding[-(64:65),] # final two rows are NAs so removing them

# read in day 1 data
day_1 <- read.csv("atsb_working_code/DB monthly natural sugar and ASB feeding/day 1-Table 1.csv")


sugar_feeding |>
  filter(Sample.female.Anth..Day.2 > 5) -> sugar_feeding
day_1 |>
  filter(Sample.female.Anth..Day.1 > 5) -> day_1
plot(sugar_feeding$females.Anth..positive.Day.2/sugar_feeding$Sample.female.Anth..Day.2,
     sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2,
     pch=20, frame.plot = F, ylab="Proportion ASB fed", xlab="Proportion anthrone positive",
     xlim=c(0,1), ylim=c(0,1))
lines(c(0,1), c(0,1))
y <- sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2
x <- sugar_feeding$females.Anth..positive.Day.2/sugar_feeding$Sample.female.Anth..Day.2
abline(lm(y ~ x), col="dodgerblue", lwd=2)

sugar_feeding |>
  mutate(anthrone_positive = females.Anth..positive.Day.2/Sample.female.Anth..Day.2,
         ASB_positive = females.ASB.positive/TOTAL.Sample.females.Day.2) -> sugar_feeding

day_1 |>
  mutate(anthrone_positive = females.Anth..positive.Day.1/Sample.female.Anth..Day.1) -> day_1

plot(sugar_feeding$anthrone_positive,
     sugar_feeding$ASB_positive,
     pch=20, frame.plot = F, ylab="Proportion ASB fed", xlab="Proportion anthrone positive",
     xlim=c(0,1), ylim=c(0,1))
lines(c(0,1), c(0,1))


# bootstrap_CI <- function(data) {
  d_anthrone <- matrix(0, nrow=9, ncol=3)
  d_asb <- matrix(0, nrow=9, ncol=3)
  for (i in 4:12) {
    sugar_feeding |>
      filter(month == i) -> t
    x <- matrix(rep(t$ASB_positive, 10000), byrow = T, ncol = length(t$ASB_positive))
    x2 <- matrix(rep(t$anthrone_positive, 10000), byrow = T, ncol = length(t$anthrone_positive))
    a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,7,T))-mean(x)})
    a2 <- apply(x2, MARGIN = 1, FUN = function(x2){mean(sample(x2,7,T))-mean(x2)})
    d_asb[i-3,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$ASB_positive)
    d_anthrone[i-3,] <- quantile(a2, c(0.025,0.5,0.975)) + mean(t$anthrone_positive)
  }
# }
  d_anthrone_day1 <- matrix(0, nrow=9, ncol=3)
  for (i in 4:12) {
    day_1 |>
      filter(month == i) -> t    
    x2 <- matrix(rep(t$anthrone_positive, 10000), byrow = T, ncol = length(t$anthrone_positive))
    a2 <- apply(x2, MARGIN = 1, FUN = function(x2){mean(sample(x2,7,T))-mean(x2)})
    d_anthrone_day1[i-3,] <- quantile(a2, c(0.025,0.5,0.975)) + mean(t$anthrone_positive)
  }

  plot(sugar_feeding$month, 
       sugar_feeding$anthrone_positive,
       pch=20, frame.plot = F, xlab="Month", ylab="Proportion positive", ylim=c(0,1))
  points(sugar_feeding$month,
         sugar_feeding$ASB_positive,
         col=2,
         pch=20)
  points(day_1$month,
         day_1$anthrone_positive,
         col=4,
         pch=20)
  lines(4:12, d_anthrone[,2], lwd=2)
  lines(4:12, d_asb[,2], lwd=2, col=2)
  lines(4:12, d_anthrone_day1[,2], lwd=2, col=4)
  polygon(c(4:12,12:4), c(d_anthrone[,1], rev(d_anthrone[,3])), 
          col=adjustcolor("black", alpha.f = 0.2), border = F)  
  polygon(c(4:12,12:4), c(d_asb[,1], rev(d_asb[,3])), 
          col=adjustcolor("red", alpha.f = 0.2), border = F)   
  polygon(c(4:12,12:4), c(d_anthrone_day1[,1], rev(d_anthrone_day1[,3])), 
          col=adjustcolor(4, alpha.f = 0.2), border = F)  
  legend(x="topright", legend=c("ASB Day 2","Anthrone Day 1", "Anthrone Day 2"), col=c(2,4,1), lty=1, lwd=2, bty="n")

  
  
  cor(d_anthrone[,2], d_asb[,2])
  
  plot(d_anthrone[,2], d_asb[,2], pch=20, frame.plot = F)
  