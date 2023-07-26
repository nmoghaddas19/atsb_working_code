insec_resis <- read.csv("~/Documents/GitHub/atsb_working_code/220415 Current insecticide resistance database for 2021_2022 ATSB study areas in Western Zambia_JC/Sheet 1-Table 1.csv")
library(dplyr)
#library(boot)
insec_resis |>
  filter(Insecticide.Class == "Pyrethroid" &
           Total.female.mosquitoes.in.all.test.replicates > 0) -> pyrethroid_resis
pyrethroid_resis$female_dead_fraction <- as.double(pyrethroid_resis$Female_Recorded.average.mortality.in.treatments....or.number.dead.)
pyrethroid_resis$female_dead_fraction[1] <- 100
pyrethroid_resis$total_females <- pyrethroid_resis$Total.female.mosquitoes.in.all.test.replicates
weighted.mean(pyrethroid_resis$female_dead_fraction, pyrethroid_resis$total_females)

table(pyrethroid_resis$Species.used.in.test.replicates)

# by species

colnames(pyrethroid_resis)[8] <- "Species"
colnames(pyrethroid_resis)
table(pyrethroid_resis$Species)
pyrethroid_resis$Species <-  ifelse(pyrethroid_resis$Species == "An. gambiae s.l", "An.gambiae s.l",
                               ifelse(pyrethroid_resis$Species == "An.funestus", "An.funestus s.l", 
                                      pyrethroid_resis$Species))
table(pyrethroid_resis$Species)

pyrethroid_resis |>
  group_by(Species) |>
  summarise(bioassay_mortality = weighted.mean(female_dead_fraction, total_females),
            error = sd(female_dead_fraction)/sqrt(n())) -> 
  species_resistance

bootstrap <- function(data, species) {
  statistic <- c()
  data |>
    filter(Species == species) -> pyrethroid_resistance
  size = length(pyrethroid_resistance$female_dead_fraction)
  for (i in 1:5000) {
    boot_sample <- sample(pyrethroid_resistance$female_dead_fraction, 
                          prob = pyrethroid_resistance$total_females, size = size,
                          replace = TRUE)
    statistic[i] <- mean(boot_sample)
  }
  return(
    list(
      statistic = statistic,
      quantiles = quantile(statistic, probs = c(0.025,0.975)),
      size = size,
      species = species,
      mean = mean(statistic)
    )
  )
}
funestus_resistance <- bootstrap(pyrethroid_resis, "An.funestus s.l")
gambiae_resistance <- bootstrap(pyrethroid_resis, "An.gambiae s.l")

