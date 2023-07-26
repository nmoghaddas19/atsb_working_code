remove.packages("malariasimulation")

p <- get_parameters(
  overrides = list(
    atsb=TRUE,
    mu_atsb = c(0.09,0.09,0.09)
  )
)
p <- set_bednets(
  p,
  timesteps = matrix(c(.533), nrow = 1, ncol = 1),
  coverages = c(.5),  # Each round is distributed to 50% of the population.
  retention = 5 * 365, # Nets are kept on average 5 years
  dn0 = matrix(c(.533), nrow = 1, ncol = 1), # Matrix of death probabilities for each mosquito species over time
  rn = matrix(c(.56), nrow = 1, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
  rnm = matrix(c(.24), nrow = 1, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
  gamman = rep(2.64 * 365, 1) # Vector of bed net half-lives for each distribution timestep
)
p <- set_atsb(parameters = p, timesteps = 50:100, coverages = rep(1,51))
out <- run_simulation(timesteps = 365, parameters = p)
