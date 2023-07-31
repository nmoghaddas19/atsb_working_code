
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

packages <- c("malariasimulation")
sources <- "ghana_noIRS_source.R"
src <- conan::conan_sources(c("nmoghaddas19/malariasimulation"))
ctx <- context::context_save("root", packages = packages, sources = sources, 
                             package_sources = src)

obs <- didehpc::queue_didehpc(ctx)


t <- obs$enqueue(
  run_country_irs()
)

config <- didehpc::didehpc_config(cores = 8)
obj <- didehpc::queue_didehpc(ctx, config = config)
s <- obj$enqueue(
  run_country_irs()
)

packages <- c("cali", "site", "grr")
sources <- "ghana_noIRS_source.R"
src <- conan::conan_sources(c("mrc-ide/cali", "mrc-ide/site"))
ctx <- context::context_save("root2", packages = packages, sources = sources, 
                             package_sources = src)
obz <- didehpc::queue_didehpc(ctx)
obz$install_packages("nmoghaddas19/malariasimulation")

z <- obz$enqueue(
  set_atsb()
)






