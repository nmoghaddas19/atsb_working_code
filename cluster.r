
library(foresite)
library(site)
library(malariasimulation)
library(grr)
library(cali)
drat:::add("mrc-ide")
# install.packages("didehpc")
library(didehpc)

setwd("Q:/")
options(didehpc.cluster = "fi--didemrchnb",
        didehpc.home = "Q:")
didehpc::didehpc_config()

packages <- c("cali")
sources <- "sources.R"
src <- conan::conan_sources(c("mrc-ide/cali"))
ctx <- context::context_save("contexts", packages = packages, sources = sources, 
                             package_sources = src)

obs <- didehpc::queue_didehpc(ctx)

t <- obj$enqueue(
  run_simulation(1000)
)
s <- obj$enqueue(
  3+4
)

t$wait(120)
t$result()

setwd("T:/")
obj <- didehpc::queue_didehpc(ctx)























# load packages local
pacman::p_load(tidyverse, didehpc, rio, cali, malariasimulation,
               malariaEquilibrium)

# Cluster setup
options(didehpc.username = "nmoghad1",
        didehpc.home = "X:")

getwd()

didehpc::didehpc_config()

drat:::add("mrc-ide")
malaria <- didehpc::path_mapping("Malaria", "X:", "\\\\fi--san03.dide.ic.ac.uk\\homes\\ddee", "X:")
malaria <- didehpc::path_mapping(name = "Malaria",
                                 "X:", 
                                 "\\fi--didef3.dide.ic.ac.uk\malaria")

config <- didehpc::didehpc_config(credentials = "ddee",
                                  shares = malaria,
                                  use_rrq =  FALSE,
                                  cluster = "fi--didemrchnb",
                                  parallel = FALSE,
                                  cores = 1)

packages <- c("dplyr", "stringr", "tidyr", "lubridate", "purrr", "mvw",
              "malariasimulation", "malariaEquilibrium", "cali")

sources = c("protopopoff_times.R", # dates and times for the mosha trial parameterisation
            "functions.R", # Isaac's small functions
            "get_params.R",  # set baseline parameters
            "sim_baseline_EIR.R", # calibrate the cluster
            "sim_forward.R") # run the model

src <- conan::conan_sources("mrc-ide/mvw") # using the malariaverse workshop site to get packages

ctx <- context::context_save(path = "new_contexts", packages = packages, sources = sources, package_sources = src)
obj = didehpc::queue_didehpc(ctx, config = config)

