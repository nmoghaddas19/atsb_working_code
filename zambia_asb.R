# The aim of thi script is to examine the ASB trial data from Zambia to 
# parameterise malariasimulation and compare predictions to observed values 
# for mosquito counts (are there in this dataset?)

library(lubridate)
# read in data
zambia <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")

# looking at the data it seems like there is no ATSB arm. all clusters had either 
# 2 or 3 ASBs deployed per eligible structure throughout so I cant compare predicted
# counts to anything if i generate them

# i guess we can still have a look at the feeding rate. first we need to collapse
# the data frame

zambia |> 
  mutate(Month=month(parse_date_time(Day.of.Date.Trap, "dmy"))) |>
  group_by(Month, Cluster) |>
  summarise(ASB_fraction=sum(form.DyeFed)/n(), 
            Count = n()) -> zambia

d_asb <- matrix(0, nrow=3, ncol=3)
for (i in 3:5) {
  zambia |>
    filter(Month == i) -> t
  x <- matrix(rep(t$ASB_fraction, 5000), byrow = T, ncol = length(t$ASB_fraction))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,length(t$ASB_fraction),T))-mean(x)})
  d_asb[i-2,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$ASB_fraction)
}  
plot(zambia$Month,
     zambia$ASB_fraction,
     pch=20,
     frame.plot = F)
polygon(c(3:5,5:3), c(d_asb[,1],rev(d_asb[,3])),
        border = FALSE, col = adjustcolor("black", alpha.f = 0.2))

d_asb <- matrix(0, nrow=3, ncol=3)
for (i in 3:5) {
  zambia |>
    filter(Month == i) -> t
  x <- matrix(rep(t$Count, 5000), byrow = T, ncol = length(t$Count))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,length(t$Count),T))-mean(x)})
  d_asb[i-2,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$Count)
}  
plot(zambia$Month,
     zambia$Count,
     pch=20,
     frame.plot = F)
polygon(c(3:5,5:3), c(d_asb[,1],rev(d_asb[,3])),
        border = FALSE, col = adjustcolor("black", alpha.f = 0.2))


  