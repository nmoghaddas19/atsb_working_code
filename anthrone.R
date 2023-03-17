# Looking at anthrone feeding and if it correlates to ASB feeding

# read in the sugar feeding data
sugar_feeding <- read.csv("atsb_working_code/DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
sugar_feeding <- sugar_feeding[-(64:65),] # final two rows are NAs so removing them

sugar_feeding |>
  filter(Sample.female.Anth..Day.2 > 10) -> sugar_feeding
plot(sugar_feeding$females.Anth..positive.Day.2/sugar_feeding$Sample.female.Anth..Day.2,
     sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2,
     pch=20, frame.plot = F, ylab="Proportion ASB fed", xlab="Proportion anthrone positive",
     xlim=c(0,0.5), ylim=c(0,1))
y <- sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2
x <- sugar_feeding$females.Anth..positive.Day.2/sugar_feeding$Sample.female.Anth..Day.2
abline(lm(y ~ x), col="dodgerblue", lwd=2)

plot(sugar_feeding$month, 
     sugar_feeding$females.Anth..positive.Day.2/sugar_feeding$Sample.female.Anth..Day.2,
     pch=20, frame.plot = F, xlab="Month", ylab="Anthrone positive")

