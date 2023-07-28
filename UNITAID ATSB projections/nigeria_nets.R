library(foresite)
library(site)
library(malariasimulation)
library(grr)
library(cali)

pyr_nets <- read.csv("~/Documents/GitHub/pyrethroid_only_nets.csv")
ig2_nets <- read.csv("~/Documents/GitHub/pyrethroid_pyrrole_nets.csv")
pbo_nets <- read.csv("~/Documents/GitHub/pyrethroid_pbo_nets.csv")

sites <- NGA$sites
site_index <- which(NGA$sites$name_1 == "Borno")
site <- single_site(NGA, site_index[i])
params <- site_parameters(
  interventions = site$interventions,
  demography = site$demography,
  vectors = site$vectors,
  seasonality = site$seasonality,
  eir = site$eir$eir[1],
  overrides = list(human_population = 10000)
)
control <- run_simulation(timesteps = params$timesteps + 5*365, 
                          parameters = params)
par(las = 1, mfrow = c(1,1))
plot(
  control$timestep/365+2000,
  control$n_detect_730_3649/control$n_730_3649,
  type = "l",
  lwd = 2,
  col = "darkorchid3",
  frame.plot = FALSE,
  ylim = c(0,1)
)
grid()
