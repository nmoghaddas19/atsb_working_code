library(devtools)
install_github("mrc-ide/odin")
install_github("mrc-ide/deterministic-malaria-model")
install.packages("dde")
library(odin)
library(ICDMM)

out <- run_model_example()
out$plot
View(out$dat)
