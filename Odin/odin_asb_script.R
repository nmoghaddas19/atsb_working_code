library(odin)
# library(ICDMM)

out <- run_model_example()
out$plot
View(out$dat)

# ITN only
out <- run_model(model = "odin_model_asb",
                 init_EIR = 50,
                 time = 1000,
                 asb_on = 500,
                 asb_off = 650,
                 feeding_rate = 0.15,
                 u_asb = 0.0,
                 admin2 = "Kayes",
                 dye_days = 0)

out2 <- run_model(model = "odin_model_asb",
                 init_EIR = 50,
                 time = 1000,
                 asb_on = 500,
                 asb_off = 650,
                 feeding_rate = 0.15,
                 u_asb = 0.0,
                 country = "Zambia",
                 admin2 = "Western",
                 dye_days = 0,
                 ITN_IRS_on = 470,
                 itn_cov = 0.0,
                 irs_cov = 0.3,
                 d_IRS0 = 1.5,
                 num_int = 4)
par(las=1)
plot(out$t,out$prev, main= "Prevalance", type='l', ylim = c(0, 1))
lines(out2$t, out2$prev, col = "blue")
abline(v = 365, lty = 2)

plot(out$t[4950:6750], (out$Svasb+out$Evasb+out$Ivasb)[4950:6750], main= "Stained",
     type='l', xlab = "Day", ylab = "No. stained mosquitoes", lwd=2,
     col=2, frame.plot = F)
lines(out2$t, out2$Svasb+out2$Evasb+out2$Ivasb, col = 1, lwd=2)
grid()
legend(x="topright", legend = c("15% bait feeding", "  0% bait feeding"),
       col=c(2, 1), lty=1, lwd=2, bty="n")

plot(out$t[4950:6750], (out$Svasb+out$Evasb+out$Ivasb)[4950:6750]/out$mv[4950:6750],
     type='l', col=1, lwd=2, frame.plot = F, xlab = "Day", ylab = "Proportion dye fed",
     ylim=c(0,1))
grid()
lines(out2$t[4950:6750], (out2$Svasb+out2$Evasb+out2$Ivasb)[4950:6750]/out2$mv[4950:6750],
      type='l', col=2, lwd=2)
legend(x="topleft", legend = c("dye longevity infinite", "dye longevity 10 days",
                               "dye longevity 4.5 days", "dye longevity 2.5 days"),
       col=c(1, 2, 4, 7), lty=1, lwd=2, bty="n")
((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5500:5550]
max((out$Svasb+out$Evasb+out$Ivasb)/out$mv)

plot(out2$t[4650:8650], out2$mu[4650:8650],
     type='l', col=2, lwd=2, frame.plot = F, xlab = "Day", ylab = "mu",
     ylim=c(0,1))

seasonality_data <- ICDMM::load_file("admin_units_seasonal.rds")

# vary feeding rate and measure dyed fraction
dyed_fraction <- c()
for (i in 1:50) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = i/100,
                   u_asb = 0.25,
                   dye_days = 1/2.5,
                   mu0 = 0.096)
  dyed_fraction[i] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5300])
}
plot(1:50/100, dyed_fraction, type="l", lwd=2, frame.plot = F,
     xlab="Feeding rate", ylab="Dyed fraction", ylim=c(0,1))
grid()

dyed_fraction_24hr <- c()
for (i in 1:50) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = i/100,
                   u_asb = 0.10)
  dyed_fraction_24hr[i] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5300])
}
lines(1:50/100, dyed_fraction_24hr, lwd=2,
     col = "dodgerblue")
legend(x="topleft", legend = c("0% mortality due to dye", "5% mortality due to dye", "10% mortality due to dye"),
       col=c(1, 2, "dodgerblue"), lty=1, lwd=2, bty="n")


# vary time of day traps are active and measure dyed fraction
dyed_fraction <- c()
for (i in 1:20) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = 0.15,
                   u_asb = 0.0)
  dyed_fraction[i] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[(4000+i):(4005+i)])
}
plot(0:19/10*24, dyed_fraction, type = "l", lwd =2, frame.plot = F, ylab="Dyed fraction",
     xlab="Hours after ASB deployment", xlim=c(0,50), ylim=c(0,max(dyed_fraction)+0.02))
grid()

dyed_fraction_10 <- c()
for (i in 1:20) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = 0.15,
                   u_asb = 0.05)
  dyed_fraction_10[i] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[(4000+i):(4005+i)])
}
lines(0:19/10*24, dyed_fraction_10, lwd =2, col="dodgerblue")
legend(x="topleft", legend = c("15% bait feeding", "10% bait feeding"),
       col=c(1, "dodgerblue"), lty=1, lwd=2, bty="n")

# vary background mortality/ ITN use
dyed_fraction <- c()
for (i in seq(0,100,5)) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = 0.10,
                   u_asb = 0.00,
                   ITN_IRS_on = 365,
                   itn_cov = i/100,
                   num_int = 2)
  dyed_fraction[i/5+1] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5300])
}
plot(seq(0,100,5), dyed_fraction, type = "l", lwd =2, frame.plot = F, ylab="Dyed fraction",
     xlab="ITN coverage", ylim=c(0,1), xlim=c(0,100))
grid()

dyed_fraction_35 <- c()
for (i in seq(0,100,5)) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = 0.05,
                   u_asb = 0.0,
                   ITN_IRS_on = 365,
                   itn_cov = i/100,
                   num_int = 2)
  dyed_fraction_35[i/5 +1] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5300])
}
lines(seq(0,100,5), dyed_fraction_35, lwd=2, col=2)
legend(x="topleft", legend = c("10% bait feeding", " 5% bait feeding"),
       col=c(1,2), lty=1, lwd=2, bty="n")
plot(out$t, out$mv, frame.plot = F, lwd=2, type="l")
plot(out$t, out$betaa, frame.plot = F, lwd=2, type="l")

# okay what about background mortality
dyed_fraction <- c()
for (i in seq(1,50,3)) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = 0.10,
                   u_asb = 0.0,
                   mu0 = i/100)
  dyed_fraction[(i+2)/3] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5300])
}
plot(seq(1,50,3), dyed_fraction, type = "l", lwd =2, frame.plot = F,
     ylab="Dyed fraction", xlab="Background mortality", ylim=c(0,1))
grid()
plot(out$t, out$mv, type = "l", lwd=2)

# varying time of year

dyed_fraction <- c()
for (i in seq(1,365,20)) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500+i,
                   asb_off = 650+i,
                   feeding_rate = 0.10,
                   u_asb = 0.0,
                   dye_days = 0)
  dyed_fraction[(i+19)/20] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5300+i*10])
}
plot(seq(1,365,20), dyed_fraction, type = "l", lwd =2, frame.plot = F,
     ylab="Dyed fraction", xlab="Offset (days)", ylim = c(0,1))
grid()
lines(out$t[5000:8650]-500, out$mv[5300:8950]/70,
     type='l', col=adjustcolor(col=2, alpha.f = 0.3), lwd=2, frame.plot = F, xlab = "Day", ylab = "No. mosquitoes")

# varying dye half life

dyed_fraction <- c()
for (i in seq(0.5,10,0.5)) {
  out <- run_model(model = "odin_model_asb",
                   init_EIR = 50,
                   time = 1000,
                   admin2 = "Kayes",
                   asb_on = 500,
                   asb_off = 650,
                   feeding_rate = 0.15,
                   u_asb = 0.0,
                   dye_days = 1/i)
  dyed_fraction[(i*2)] <- mean(((out$Svasb+out$Evasb+out$Ivasb)/out$mv)[5300])
}
plot(seq(0.5,10,0.5), dyed_fraction, type = "l", lwd =2, frame.plot = F,
     ylab="Dyed fraction", xlab="Dye longevity (days)", ylim=c(0,0.5), xlim=c(0,10))
grid()
plot(out$t, out$mv, type = "l", lwd=2)
