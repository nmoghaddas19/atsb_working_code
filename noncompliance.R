################################################################################
## The purpose of this script is to characterise what is the epidemiological  ##
## effect of having two mosquito populations that differ only in their        ##
## propensity to feed on sugar baits. For example what is the difference      ##
## between have two populations which each have a feeding rate of 10% vs      ##
## having one with 5% and another with 15%?                                   ##
################################################################################
detach("package:malariasimulation", unload = TRUE)
library(malariasimulation)
library(site)
library(foresite)
library(RColorBrewer)
library(DescTools)
################################################################################
############################### Scenario 1: 10%  ###############################
################################################################################

mali <- MLI
kayes_rural <- single_site(mali, 5)
feed <- c(0.10, 0.15, 0.20)
# props <- c(1, 0.75, 0.5)
out_nc <- list()
for (i in 1:length(feed)) {
  name <- as.character(feed[i])
  kayes_rural_params <- site_parameters(
    interventions = kayes_rural$interventions,
    demography = kayes_rural$demography,
    vectors = kayes_rural$vectors,
    seasonality = kayes_rural$seasonality,
    eir = kayes_rural$eir$eir[1],
    overrides = list(human_population = 5000,
                     mu_atsb = c(feed[i],0.2-feed[i]))
  )
  arab_params_1 <- gamb_params
  arab_params_1$species <- "arab"
  kayes_rural_params <- set_species(kayes_rural_params, species=list(gamb_params,arab_params_1),
                                    proportions=c(0.5,0.5))
  kayes_rural_params <- set_atsb(parameters = kayes_rural_params,
                                 timesteps = (17*365+5*30):(18*365), 
                                 coverages = rep(1,366-5*30))
  out_nc[[name]] <- run_simulation(timesteps = kayes_rural_params$timesteps,
                                       parameters = kayes_rural_params)
  print(i)
}

par(las=1)
plot(out_nc[[1]]$timestep/365+2000,
     out_nc[[1]]$total_M_gamb,
     type="l", lwd=2, frame.plot = F, xlab="Year", ylab="Population", xlim=c(2016,2018))
lines(out_nc[[1]]$timestep/365+2000,
     out_nc[[1]]$total_M_arab,
     lty=2, lwd=2)
lines(out_nc[[2]]$timestep/365+2000,
      out_nc[[2]]$total_M_gamb,
      col=2, lwd=2)
lines(out_nc[[2]]$timestep/365+2000,
      out_nc[[2]]$total_M_arab,
      col=2, lwd=2)
legend(x="topright", legend=c("10% and 10%", "5% and 15%"), col=c(1,2), 
       lty = 1, lwd=2, bty = "n")


plot(out_nc[[1]]$timestep/365 +2000, 
     out_nc[[1]]$total_M_arab + out_nc[[1]]$total_M_gamb,
     type="l", lwd=2, frame.plot = F, xlim=c(2016,2018), xlab="Year", ylab=NA,
     cex.axis=1.2)
lines(out_nc[[2]]$timestep/365+2000, 
      out_nc[[2]]$total_M_arab + out_nc[[2]]$total_M_gamb,
      lwd=2, col=2)
lines(out_nc[[3]]$timestep/365+2000, 
      out_nc[[3]]$total_M_arab + out_nc[[3]]$total_M_gamb,
      lwd=2, col=4)
legend(x="topright", legend=c("10% and 10%", "5% and 15%", "0% and 20%"), col=c(1,2,4), 
       lty = 1, lwd=2, bty = "n", cex=1.2)

plot(out_nc[[4]]$timestep/365+2000, 
     out_nc[[4]]$n_detect_730_3649 / out_nc[[4]]$n_730_3649, ylim = c(0,1),
     type="l", lwd=2, frame.plot = F, xlim=c(2016,2020), xlab="Year", ylab=NA, 
     cex.axis=1.2)
lines(out_nc[[1]]$timestep/365+2000, 
      out_nc[[1]]$n_detect_730_3649 / out_nc[[1]]$n_730_3649,
      lwd=2, col=7)
lines(out_nc[[2]]$timestep/365+2000, 
      out_nc[[2]]$n_detect_730_3649 / out_nc[[2]]$n_730_3649,
      lwd=2, col=2)
lines(out_nc[[3]]$timestep/365+2000, 
      out_nc[[3]]$n_detect_730_3649 / out_nc[[3]]$n_730_3649,
      lwd=2, col=4)
legend(x="topright", legend=c("10% and 10%", "5% and 15%", "0% and 20%", "No ATSB"), col=c(7,2,4,1), 
       lty = 1, lwd=2, bty = "n", cex=1.2)

kayes_rural_params <- site_parameters(
  interventions = kayes_rural$interventions,
  demography = kayes_rural$demography,
  vectors = kayes_rural$vectors,
  seasonality = kayes_rural$seasonality,
  eir = kayes_rural$eir$eir[1],
  overrides = list(human_population = 5000)
)
arab_params_1 <- gamb_params
arab_params_1$species <- "arab"
kayes_rural_params <- set_species(kayes_rural_params, species=list(gamb_params,arab_params_1),
                                  proportions=c(0.5,0.5))
out_nc[["0"]] <- run_simulation(timesteps = kayes_rural_params$timesteps,
                                 parameters = kayes_rural_params)

################################################################################
############################### Scenario 2: 35%  ###############################
################################################################################

mali <- MLI
kayes_rural <- single_site(mali, 5)
feed <- c(0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.7)
out_nc_35 <- list()
for (i in 1:length(feed)) {
  name <- as.character(feed[i])
  kayes_rural_params <- site_parameters(
    interventions = kayes_rural$interventions,
    demography = kayes_rural$demography,
    vectors = kayes_rural$vectors,
    seasonality = kayes_rural$seasonality,
    eir = kayes_rural$eir$eir[1],
    overrides = list(human_population = 5000,
                     mu_atsb = c(feed[i],0.70-feed[i]))
  )
  arab_params_1 <- gamb_params
  arab_params_1$species <- "arab"
  kayes_rural_params <- set_species(kayes_rural_params, species=list(gamb_params,arab_params_1),
                                    proportions=c(0.5,0.5))
  kayes_rural_params <- set_atsb(parameters = kayes_rural_params,
                                 timesteps = (17*365+5*30):(18*365), 
                                 coverages = rep(1,366-5*30))
  out_nc_35[[name]] <- run_simulation(timesteps = kayes_rural_params$timesteps,
                                   parameters = kayes_rural_params)
  print(i)
}

plot(out_nc_35[[1]]$timestep/365+2000,
     out_nc_35[[1]]$total_M_gamb,
     type="l", lwd=2, frame.plot = F, xlab="Year", ylab="Population", xlim=c(2016,2018))
lines(out_nc_35[[1]]$timestep/365+2000,
      out_nc_35[[1]]$total_M_arab,
      lty=2, lwd=2)
lines(out_nc_35[[2]]$timestep/365+2000,
      out_nc_35[[2]]$total_M_gamb,
      col=2, lwd=2)
lines(out_nc_35[[2]]$timestep/365+2000,
      out_nc_35[[2]]$total_M_arab,
      col=2, lwd=2)
legend(x="topright", legend=c("10% and 10%", "5% and 15%"), col=c(1,2), 
       lty = 1, lwd=2, bty = "n")


plot(out_nc_35[[1]]$timestep/365 +2000, 
     out_nc_35[[1]]$total_M_arab + out_nc_35[[1]]$total_M_gamb,
     type="l", lwd=2, frame.plot = F, xlim=c(2016,2018), xlab="Year", ylab=NA,
     cex.axis=1.2, col=brewer.pal(n=8,name="RdYlBu")[8])
for (i in 2:length(out_nc_35)) {
  lines(out_nc_35[[i]]$timestep/365+2000, 
        out_nc_35[[i]]$total_M_arab + out_nc_35[[i]]$total_M_gamb,
        lwd=2, col=brewer.pal(n=8,name="RdYlBu")[i])
}
legend(x="topleft", legend=rev(c("35% and 35%", "30% and 40%", "25% and 45%", "20% and 50%",
                              "15% and 55%", "10% and 60%","5% and 65%", "0% and 70%", "No ATSB")), 
       col= rev(c(brewer.pal(n=8,name="RdYlBu"),"black")), lty = 1, lwd=2, bty = "n")
lines(out_nc[[4]]$timestep/365+2000, 
      out_nc[[4]]$total_M_arab + out_nc[[4]]$total_M_gamb,
      lwd=2, col="black")

?ColorBrewer
plot(out_nc_35[[1]]$timestep/365+2000, 
     out_nc_35[[1]]$n_detect_730_3649 / out_nc_35[[1]]$n_730_3649, ylim = c(0,1),
     type="l", lwd=2, frame.plot = F, xlim=c(2016,2020), xlab="Year", ylab=NA,
     cex.axis=1.2, col=brewer.pal(n=8,name="RdYlBu")[8])
for (i in 2:length(out_nc_35)) {
  lines(out_nc_35[[i]]$timestep/365+2000, 
        out_nc_35[[i]]$n_detect_730_3649 / out_nc_35[[i]]$n_730_3649,
        lwd=2, col=brewer.pal(n=8,name="RdYlBu")[i])
}
legend(x="topright", legend=rev(c("35% and 35%", "30% and 40%", "25% and 45%", "20% and 50%",
                              "15% and 55%", "10% and 60%","5% and 65%", "0% and 70%", "No ATSB")), 
       col= rev(c(brewer.pal(n=8,name="RdYlBu"),"black")), lty = 1, lwd=2, bty = "n")
lines(out_nc[[4]]$timestep/365+2000, 
      out_nc[[4]]$n_detect_730_3649 / out_nc[[4]]$n_730_3649,
      lwd=2, col="black")

plot(out_nc_35[[1]]$timestep/365+2000, 
     out_nc_35[[1]]$n_detect_730_3649,
     type="l", lwd=2, frame.plot = F, xlim=c(2016,2020), xlab="Year", ylab=NA,
     cex.axis=1.2, col=brewer.pal(n=8,name="RdYlBu")[8])
for (i in 2:length(out_nc_35)) {
  lines(out_nc_35[[i]]$timestep/365+2000, 
        out_nc_35[[i]]$n_detect_730_3649,
        lwd=2, col=brewer.pal(n=8,name="RdYlBu")[i])
}
legend(x="topright", legend=rev(c("35% and 35%", "30% and 40%", "25% and 45%", "20% and 50%",
                                  "15% and 55%", "10% and 60%","5% and 65%", "0% and 70%", "No ATSB")), 
       col= rev(c(brewer.pal(n=8,name="RdYlBu"),"black")), lty = 1, lwd=2, bty = "n")
lines(out_nc[[4]]$timestep/365+2000, 
      out_nc[[4]]$n_detect_730_3649,
      lwd=2, col="black")

cases_averted <- c()
for (i in 1:length(out_nc_35)) {
  cases_averted[i] <- AUC( (17*365+5*30):(20*365+5*30), -out_nc_35[[i]]$n_detect_730_3649[(17*365+5*30):(20*365+5*30)]
                            +(out_nc[[4]]$n_detect_730_3649)[(17*365+5*30):(20*365+5*30)])
}

par(las=1)
plot(NA)
barplot(rev(c(cases_averted/3)), ylim=c(0,20000),
        ylab=NA,
        col=brewer.pal(n=8,name="RdYlBu"))

text(x = 1:9,
     labels = rev(c("35% and 35%", "30% and 40%", "25% and 45%", "20% and 50%",
                              "15% and 55%", "10% and 60%","5% and 65%", "0% and 70%")),
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 1)
