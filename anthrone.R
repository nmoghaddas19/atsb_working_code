# Looking at anthrone feeding and if it correlates to ASB feeding

library(dplyr)
library(boot)

# read in the sugar feeding data
sugar_feeding <- read.csv("atsb_working_code/DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
sugar_feeding <- sugar_feeding[-(64:65),] # final two rows are NAs so removing them

sugar_feeding |>
  filter(Sample.female.Anth..Day.2 > 10) -> sugar_feeding
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

plot(sugar_feeding$anthrone_positive,
     sugar_feeding$ASB_positive,
     pch=20, frame.plot = F, ylab="Proportion ASB fed", xlab="Proportion anthrone positive",
     xlim=c(0,1), ylim=c(0,1))
lines(c(0,1), c(0,1))


# bootstrap_CI <- function(data) {
  d_anthrone <- matrix(0, nrow=7, ncol=3)
  d_asb <- matrix(0, nrow=7, ncol=3)
  for (i in 6:12) {
    sugar_feeding |>
      filter(month == i) -> t
    x <- matrix(rep(t$ASB_positive, 10000), byrow = T, ncol = length(t$ASB_positive))
    x2 <- matrix(rep(t$anthrone_positive, 10000), byrow = T, ncol = length(t$anthrone_positive))
    a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,7,T))-mean(x)})
    a2 <- apply(x2, MARGIN = 1, FUN = function(x2){mean(sample(x2,7,T))-mean(x2)})
    d_asb[i-5,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$ASB_positive)
    d_anthrone[i-5,] <- quantile(a2, c(0.025,0.5,0.975)) + mean(t$anthrone_positive)
  }
# }

  plot(sugar_feeding$month, 
       sugar_feeding$anthrone_positive,
       pch=20, frame.plot = F, xlab="Month", ylab="Proportion positive", ylim=c(0,1))
  points(sugar_feeding$month,
         sugar_feeding$ASB_positive,
         col=2,
         pch=20)
  lines(6:12, d_anthrone[,2], lwd=2)
  lines(6:12, d_asb[,2], lwd=2, col=2)
  polygon(c(6:12,12:6), c(d_anthrone[,1], rev(d_anthrone[,3])), 
          col=adjustcolor("black", alpha.f = 0.2), border = F)  
  polygon(c(6:12,12:6), c(d_asb[,1], rev(d_asb[,3])), 
          col=adjustcolor("red", alpha.f = 0.2), border = F)   
  legend(x="topright", legend=c("ASB","Anthrone"), col=c(2,1), lty=1, lwd=2, bty="n")
  
  cor(d_anthrone[,2], d_asb[,2])
  
  plot(d_anthrone[,2], d_asb[,2], pch=20, frame.plot = F)
  