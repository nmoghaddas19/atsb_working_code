library(umbrella)

# read in data
zambia <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")
names <- unique(zambia$Cluster)
zambia |> 
  mutate(Date=parse_date_time(Day.of.Date.Trap, "dmy")) |>
  group_by(Date, Cluster) |>
  summarise(ASB_fraction=sum(form.DyeFed)/n(), 
            Count = n()) -> zambia
zambia$days <- interval(parse_date_time("01-01-2021", "dmy"), zambia$Date) |> 
  as.numeric('days')
par(mfrow=c(2,5))
for (i in 1:length(names)) {
  zambia |>
    filter(Cluster == names[i]) -> this_zambia
  plot(this_zambia$days, this_zambia$Count, frame.plot = F, pch = 20, xlim=c(1,365))
  
  fit1 <- fit_fourier(rainfall = this_zambia$Count, t = this_zambia$days, floor = 0)
  predict_1 <- fourier_predict(coef = fit1$coefficients, t = 1:365, floor = 0)
  lines(predict_1, col = "deeppink", lwd = 2)
  title(main=names[i])
}






fit2 <- fit_fourier(rainfall = zambia$Count, t = zambia$days, floor = 0.001)
predict_2 <- fourier_predict(coef = fit2$coefficients, t = 1:365, floor = 0.001)
lines(predict_2, col = "dodgerblue", lwd = 2, lty = 2)
