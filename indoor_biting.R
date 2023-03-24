# test the effect of changing indoor biting
phi <- c(0.9,0.8,0.7,0.6,0.5)
mali <- MLI
kayes_rural <- single_site(mali, 5)
out_phi <- list()
for (i in 1:length(phi)) {
  name <- as.character(phi[i])
  kayes_rural_params <- site_parameters(
    interventions = kayes_rural$interventions,
    demography = kayes_rural$demography,
    vectors = kayes_rural$vectors,
    seasonality = kayes_rural$seasonality,
    eir = kayes_rural$eir$eir[1],
    overrides = list(human_population = 5000,
                     phi_indoors = c(phi[i],phi[i],phi[i]))
  )
  out_phi[[name]] <- run_simulation(timesteps = kayes_rural_params$timesteps,
                                       parameters = kayes_rural_params)
  print(i)
}

plot(out_phi[[1]]$timestep/365 +2000, 
     out_phi[[1]]$total_M_arab + out_phi[[1]]$total_M_gamb + out_phi[[1]]$total_M_fun,
     type="l", lwd=2, frame.plot = F, xlim=c(2016,2018), xlab="Year", ylab="Population")
for (i in 2:length(out_phi)) {
  lines(out_phi[[i]]$timestep/365+2000, 
        out_phi[[i]]$total_M_arab + out_phi[[i]]$total_M_gamb + out_phi[[i]]$total_M_fun,
        lwd=2, col=i)
}
legend(x="topright", legend=names(out_phi), 
       col= 1:length(out_phi), lty = 1, lwd=2, bty = "n")


plot(out_phi[[1]]$timestep/365+2000, 
     out_phi[[1]]$n_detect_730_3649 / out_phi[[1]]$n_730_3649, ylim = c(0,1),
     type="l", lwd=2, frame.plot = F, xlim=c(2016,2020), xlab="Year", ylab="PfPr2-10")
for (i in 2:length(out_phi)) {
  lines(out_phi[[i]]$timestep/365+2000, 
        out_phi[[i]]$n_detect_730_3649 / out_phi[[i]]$n_730_3649,
        lwd=2, col=i)
}
legend(x="topright", legend=names(out_phi), 
       col= 1:length(out_phi), lty = 1, lwd=2, bty = "n")
lines(out_nc[[4]]$timestep/365+2000, 
      out_nc[[4]]$n_detect_730_3649 / out_nc[[4]]$n_730_3649,
      lwd=2, col="darkgreen")