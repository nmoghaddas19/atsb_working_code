



out <- run_model(model = "odin_model_asb",
                 init_EIR = 50,
                 time = 1000,
                 asb_on = 345+60,
                 asb_off = 345+60+90,
                 feeding_rate = 0.10,
                 u_asb = 0.0,
                 country = "Zambia",
                 admin2 = "Western",
                 dye_days = 1/7,
                 ITN_IRS_on = 365,
                 itn_cov = 0.6,
                 irs_cov = 0.25,
                 num_int = 4)
plot(out$t,out$prev, type='l', ylim = c(0, 1), frame.plot = F, xlab="Day", ylab="Prevalence")
lines(out$t,out$mv/60, type='l', col=adjustcolor(col=2, alpha.f = 0.5))
abline(v=365+61, lty=2)
polygon(c(365+60,365+60+90, 365+60+90, 365+60), c(1,1,0,0),
        col=adjustcolor(col=4,alpha.f = 0.4), border = F)

plot(out$t[3000:6500], (out$Svasb+out$Evasb+out$Ivasb)[3000:6500]/out$mv[3000:6500],
     type='l', col=1, lwd=2, frame.plot = F, xlab = "Day", ylab = "Proportion dye fed",
     ylim=c(0,1), xlim=c(300,650))
grid()

dyed_fraction <- matrix(0, nrow=1000, ncol=30)
plot(NULL, xlim=c(0,0.5), ylim=c(0,1), xlab="Feeding rate", ylab="Dyed fraction",
     frame.plot=F)
grid()
for (j in 101:200) {
  ITN_coverage <- runif(1, 0.4, 0.8)
  IRS_coverage <- 0
  dye_days <- runif(1, 4, 12)
  offset <- runif(1, 20, 90)
  dye_mortality <- runif(1, 0, 0.1)
  for (i in seq(1,50,2)) {
    out <- suppressMessages(run_model(model = "odin_model_asb",
                     init_EIR = 50,
                     time = 471,
                     country = "Zambia",
                     admin2 = "Western",
                     asb_on = 320,
                     asb_off = 410,
                     feeding_rate = i/100,
                     u_asb = dye_mortality,
                     dye_days = 1/dye_days,
                     ITN_IRS_on = 300,
                     itn_cov = ITN_coverage,
                     irs_cov = IRS_coverage,
                     num_int = 4))
    dyed_fraction[j,(i+5)/2+3] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[3200+offset*10])
  }
  lines(seq(1,50,2)/100, dyed_fraction[j,6:30], type="l", lwd=1.5, col=adjustcolor(col=4, alpha.f = 0.5))
  dyed_fraction[j,1] <- ITN_coverage
  dyed_fraction[j,2] <- IRS_coverage
  dyed_fraction[j,3] <- dye_days
  dyed_fraction[j,4] <- offset
  dyed_fraction[j,5] <- dye_mortality
  print(j)
}

plot(1:50/100, dyed_fraction, type="l", lwd=2, frame.plot = F,
     xlab="Feeding rate", ylab="Dyed fraction", ylim=c(0,1))
grid()

write.csv(x = dyed_fraction, file = "~/Documents/GitHub/atsb_working_code/dyed_fraction.csv")

# whats the minimum/maximum I would expect?

plot(NULL, xlim=c(0,0.5), ylim=c(0,1), xlab="Feeding rate", ylab="Dyed fraction",
     frame.plot=F)
grid()

ITN_coverage <- 0.8
IRS_coverage <- 0.5
dye_days <- 3
offset <- 20
dye_mortality <- 0.1

dyed_fraction <- numeric()
for (i in seq(1,50,2)) {
  out <- suppressMessages(run_model(model = "odin_model_asb",
                                    init_EIR = 50,
                                    time = 471,
                                    country = "Zambia",
                                    admin2 = "Western",
                                    asb_on = 320,
                                    asb_off = 470,
                                    feeding_rate = i/100,
                                    u_asb = dye_mortality,
                                    dye_days = 1/dye_days,
                                    ITN_IRS_on = 300,
                                    itn_cov = ITN_coverage,
                                    irs_cov = IRS_coverage,
                                    num_int = 4))
  dyed_fraction[(i+1)/2] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[3200+offset*10])
}
lines(seq(1,50,2)/100, dyed_fraction, type="l", lwd=1.5,
      col=adjustcolor(col=4, alpha.f = 0.5))
dyed_fraction_min <- dyed_fraction


plot(NULL, xlim=c(0,0.5), ylim=c(0,1), xlab="Feeding rate", ylab="Dyed fraction",
     frame.plot=F)
grid()
polygon(c(seq(1,50,2)/100, rev(seq(1,50,2)/100)),
        c(dyed_fraction_min,rev(dyed_fraction_max)),
        col = adjustcolor(col = 4, alpha.f = 0.5), border = F)
