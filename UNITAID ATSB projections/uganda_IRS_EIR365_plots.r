country_results <- readRDS("~/Documents/uganda_IRS_EIR365.RDS")

par(las = 1, mfrow = c(3,3))
for (i in 1:length(country_results)) {
  plot(
    country_results[[i]]$bells_whistles_0.20$timestep/365 + 2000,
    country_results[[i]]$bells_whistles_0.20$n_detect_730_3649/
      country_results[[i]]$bells_whistles_0.20$n_730_3649,
    type = "l",
    lwd = 2,
    frame.plot = FALSE,
    col = "darkorchid4",
    xlab = "Year",
    ylab = "Prevalence",
    xlim = c(2022, 2025),
    ylim = c(0, 1)
  )
  grid()
  polygon(x = c(country_results[[i]]$bells_whistles_0.20$timestep/365 + 2000,
                rev(country_results[[i]]$bells_whistles_0.20$timestep/365 + 2000)),
          y = c(country_results[[i]]$bells_whistles_0.05$n_detect_730_3649/
                  country_results[[i]]$bells_whistles_0.05$n_730_3649,
                rev(country_results[[i]]$bells_whistles_0.35$n_detect_730_3649/
                      country_results[[i]]$bells_whistles_0.35$n_730_3649)),
          border = FALSE,
          col = adjustcolor(col = "darkorchid4", alpha.f = 0.3))
  
  lines(
    country_results[[i]]$atsb_stopirs_0.20$timestep/365 + 2000,
    country_results[[i]]$atsb_stopirs_0.20$n_detect_730_3649/
      country_results[[i]]$atsb_stopirs_0.20$n_730_3649,
    lwd = 2,
    col = "mediumseagreen"
  )
  
  polygon(x = c(country_results[[i]]$atsb_stopirs_0.20$timestep/365 + 2000,
                rev(country_results[[i]]$atsb_stopirs_0.20$timestep/365 + 2000)),
          y = c(country_results[[i]]$atsb_stopirs_0.05$n_detect_730_3649/
                  country_results[[i]]$atsb_stopirs_0.05$n_730_3649,
                rev(country_results[[i]]$atsb_stopirs_0.35$n_detect_730_3649/
                      country_results[[i]]$atsb_stopirs_0.35$n_730_3649)),
          border = FALSE,
          col = adjustcolor(col = "mediumseagreen", alpha.f = 0.3))
  
  lines(
    country_results[[i]]$control$timestep/365 + 2000,
    country_results[[i]]$control$n_detect_730_3649/
      country_results[[i]]$control$n_730_3649,
    lwd = 2,
    col = "black"
  )
  
  # legend("topleft", legend = c("Continue IRS", "- IRS + ATSB + Nets", "Continue IRS + ATSB + IG2"), 
  #        col = c("black", "mediumseagreen", "darkorchid4"), lty = 1, bty = "n", 
  #        lwd = 2)
  title(country_results[[i]]$site)
}



