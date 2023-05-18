################################################################################
# In this script, I have written a workflow that takes a dataframe of the      #
# sort we are expecting from the ATSB trials and generates predictions with    #
# confidence intervals based on malariasimulation. The first step is to        #
# generate estimates of dyed fraction and uncertainty from the data. The       #
# second step is to convert the dyed fraction to daily feeding rates with.     #
# uncertainty. The third step is to plug those bounds into malariasim to       #
# generate predictions.                                                        #
################################################################################

# load in packages
library(malariasimulation)
library(ICDMM)
library(site)
library(foresite)
library(lubridate)
library(RColorBrewer)
library(lme4)

# functions
function(p, V, method = "mcculloch", inverse = FALSE) {
  if(method == "mcculloch" & inverse) {
    method <- "zeger"
    warning("The McCulloch method can't be used when inverse = TRUE. Changing to Zeger.")
  }
  stopifnot(!(method == "mcculloch" & inverse))
  Beta <- qlogis(p)
  if(method == "mcculloch") {
    return(plogis(Beta - 0.5 * V * tanh(Beta * (1 + 2 * exp(-0.5 * V))/6)))
  }
  if(method == "zeger") {
    if (inverse) {
      plogis(Beta * sqrt(256 * V / (75 * pi^2) + 1))
    } else {
      plogis(Beta/sqrt(1 + ((16 * sqrt(3))/(15 * pi))^2 * V))
    }
  }
}

# load in data
dat.all <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")
# cleaning and checking ####
names(dat.all) <- tolower(names(dat.all))
dim(dat.all)
table(table(dat.all$form.sampleid.check))
dat.all <- dat.all[!duplicated(dat.all$form.sampleid.check), ]
table(table(dat.all$form.sampleid.check))
dat.all$mosquito_species[dat.all$an..gambiae %in% 1] <- "gambiae" # make species a variable
dat.all$mosquito_species[dat.all$an..funestus %in% 1] <- "funestus" # make species a variable
dat.all$mosquito_species <- factor(dat.all$mosquito_species) # make species a factor
dat.all$location <- factor(dat.all$indoor, 0:1, c("Outdoors", "Indoors")) # make location a factor
dat.all$anoph_sex <- factor(dat.all$male, 0:1, c("female", "male")) # make sex a factor
dat.all$n.asb <- dat.all$asbs.deployed 
table(dat.all$n.asb, exclude = NULL)
dat.all$date.id <- as.Date(dat.all$day.of.date.trap, "%d-%b-%y")
dat.all$month_caught <- factor(month(dat.all$date.id))
dat.all$collection_date <- factor(dat.all$date.id)
dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 2] <-
  paste0(dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 2], ".23")
dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 3] <-
  paste0(dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 3], ".23")
dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 3] <-
  paste0(dat.all$cluster[dat.all$phase == 1 & dat.all$n.asb == 3], ".32")
dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 2] <-
  paste0(dat.all$cluster[dat.all$phase == 2 & dat.all$n.asb == 2], ".32")

dat.all$an..gambiae <- dat.all$an..funestus <- dat.all$indoor <- dat.all$outdoor <-
  dat.all$male <- dat.all$female <- dat.all$asbs.deployed <-
  dat.all$day.of.date.trap <- dat.all$cluster.code <- NULL
sum(is.na(dat.all))
colSums(is.na(dat.all))
table(dat.all$form.bloodfed, exclude = NULL)
dat.all <-
  dat.all[dat.all$mosquito_species %in% c("gambiae", "funestus") &
            dat.all$anoph_sex %in% "female", ]
table(dat.all$cluster)
table(dat.all$mosquito_species, exclude = NULL)
table(dat.all$collection_method, exclude = NULL)
table(dat.all$anoph_sex, exclude = NULL)
range(dat.all$date.id)
tapply(dat.all$date.id, dat.all$phase, range)
table(dat.all$cluster, dat.all$n.asb)
table(dat.all$cluster, dat.all$n.asb, dat.all$phase)
dat.all <-
  dat.all[order(dat.all$n.asb, dat.all$cluster, dat.all$collection_date,
                dat.all$mosquito_species, dat.all$location), ]
dat.all$cluster <- factor(dat.all$cluster, unique(dat.all$cluster))
dat.all$n.asb <- factor(dat.all$n.asb, 2:3, c("2 stations", "3 stations"))
dat.all$hh <- factor(dat.all$form.collectionid.check)
dat.all$form.collectionid.check <- NULL
table(form.bloodfed = dat.all$form.bloodfed, form.dyefed = dat.all$form.dyefed, 
      exclude = NULL)
prop.table(table(form.bloodfed = dat.all$form.bloodfed, 
                 form.dyefed = dat.all$form.dyefed, exclude = NULL))
dat.all$positive <- dat.all$form.dyefed
dat.all$form.dyefed <- NULL
dat.all$cluster.day <- factor(paste(dat.all$cluster, dat.all$date.id, sep = ":"))

# stop ####

# fit regression model for each species 
sp <- c("funestus")
dat <- droplevels(dat.all[dat.all$mosquito_species  %in% sp, ])
term <- "cluster"
rest.of.formula <- "positive ~ (1 | collection_date) + (1|hh)"
form <- paste(rest.of.formula, term, sep = " + ")                             
fit <- glmer(form, family = "binomial", data = dat, 
             control = glmerControl(optimizer = "bobyqa"))
print(summary(fit))

fit0 <- update(fit, ~ . -1)
est.tab.bias <-
  cbind(fixef(fit0),
        confint(fit0, method = "Wald")[paste0(term, levs), ])
# adjust predictions for Jensen's inequality
est.tab_fun <-
  jensen.logit.adjust(p = plogis(est.tab.bias), V = sum(unlist(VarCorr(fit0))))

# dat.all$form.bloodfed <- factor(dat.a;;)
sp <- c("gambiae")
dat <- droplevels(dat.all[dat.all$mosquito_species  %in% sp, ])
term <- "cluster"
rest.of.formula <- "positive ~ (1 | collection_date)"
form <- paste(rest.of.formula, term, sep = " + ")                             
fit <- glmer(form, family = "binomial", data = dat, 
             control = glmerControl(optimizer = "bobyqa"))
print(summary(fit))

fit0 <- update(fit, ~ . -1)
est.tab.bias <-
  cbind(fixef(fit0),
        confint(fit0, method = "Wald")[paste0(term, levs), ])
# adjust predictions for Jensen's inequality
est.tab_gamb <-
  jensen.logit.adjust(p = plogis(est.tab.bias), V = sum(unlist(VarCorr(fit0))))

dyed_fraction_bounds <- cbind(est.tab_fun, est.tab_gamb)
colnames(dyed_fraction_bounds)[c(1,4)] <- c("fun", "gamb")

# convert to daily feeding rate
dyed_fraction <- read.csv(file="~/Documents/GitHub/atsb_working_code/dyed_fraction.csv")[1:200,]
dyed_fraction <- apply(dyed_fraction, MARGIN = 2, FUN = quantile, probs = c(0.025,0.975))[,7:31]
upper_fun <- approx(x = dyed_fraction[1,], y = seq(1,50,2)/100,
                xout = dyed_fraction_bounds[,3])$y
lower_fun <- approx(x = dyed_fraction[2,], y = seq(1,50,2)/100, 
                xout = dyed_fraction_bounds[,2])$y

upper_gamb <- approx(x = dyed_fraction[1,], y = seq(1,50,2)/100,
                     xout = dyed_fraction_bounds[,6])$y
lower_gamb <- approx(x = c(0,dyed_fraction[2,]), y = c(0,seq(1,50,2)/100), 
                     xout = dyed_fraction_bounds[,5])$y

feed_rates <- cbind(control = numeric(10), lower_fun, upper_fun, lower_gamb, upper_gamb)
rownames(feed_rates) <- names(fixef(fit0))

# forecast predictions using malariasimulation
zambia <- ZMB
zambia$sites$name_1
western_rural <- single_site(zambia, 18)
# plot(western_rural$prevalence$year, western_rural$prevalence$pfpr, type="l", lwd=2, 
#      frame.plot = F, ylim = c(0,1))
# western_rural$vectors

out_data <- list()
for (j in 1:2) {
  bound <- c("lower", "upper")[j]
  for (i in 1:nrow(feed_rates)) {
    cluster <- names(fixef(fit0))[i]
    western_rural_params <- site_parameters(
      interventions = western_rural$interventions,
      demography = western_rural$demography,
      vectors = western_rural$vectors,
      seasonality = western_rural$seasonality,
      eir = western_rural$eir$eir[1],
      overrides = list(human_population = 10000,
                       mu_atsb = c(feed_rates[i,j+3], feed_rates[i,j+1], feed_rates[i,j+3]))
    )
    western_rural_params <- set_atsb(parameters = western_rural_params,
                                     timesteps = (21*365+60):(21*365+148), 
                                     coverages = rep(1,148-59))
    name <- paste(cluster, bound, sep = "_")
    out_data[[name]] <- run_simulation(timesteps = western_rural_params$timesteps,
                                          parameters = western_rural_params)
    print(i/nrow(feed_rates))
  }
  print(j)
}
western_rural_params <- site_parameters(
  interventions = western_rural$interventions,
  demography = western_rural$demography,
  vectors = western_rural$vectors,
  seasonality = western_rural$seasonality,
  eir = western_rural$eir$eir[1],
  overrides = list(human_population = 10000,
                   mu_atsb = c(0,0,0))
)
western_rural_params <- set_atsb(parameters = western_rural_params,
                                 timesteps = (21*365+60):(21*365+148), 
                                 coverages = rep(1,148-59))
for (i in 1:10) {
  if (i==1) {
    out_data[["Control"]] <- run_simulation(timesteps = western_rural_params$timesteps,
                                       parameters = western_rural_params)
  } else {
    name <- paste0("Control","_ ", i)
    out_data[[name]] <- run_simulation(timesteps = western_rural_params$timesteps,
                                       parameters = western_rural_params)
  }
  print(i/10)
}

plot(out_data$Control$timestep/365 +2000, 
     out_data$Control$total_M_arabiensis+out_data$Control$total_M_funestus+out_data$Control$total_M_gambiae,
     type="l", frame.plot = F, xlim = c(2019,2022),
     xlab = "Year", ylab = "Total mosquitoes")
# i=2
# total_M_lower <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
# total_M_upper <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
# polygon(c(out_data[[i]]$timestep/365 +2000, rev(out_data[[i+10]]$timestep/365 +2000)),
#         c(total_M_lower,rev(total_M_upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
#         border = F)
for (i in 1:10) {
  total_M_lower <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
  total_M_upper <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
  total_M_control <- out_data[[i+20]]$total_M_arabiensis + out_data[[i+20]]$total_M_funestus + out_data[[i+20]]$total_M_gambiae
  lines(out_data[[i]]$timestep/365+2000, total_M_lower, col = brewer.pal(10, "RdYlBu")[2])
  lines(out_data[[i+10]]$timestep/365+2000, total_M_upper, col = brewer.pal(10, "RdYlBu")[9])
  lines(out_data[[i+20]]$timestep/365+2000, total_M_control, col="black")
}
lower <- matrix(0,nrow = 8395, ncol=10)
upper <- matrix(0,nrow = 8395, ncol=10)
control <- matrix(0,nrow = 8395, ncol=10)
for (i in 1:10) {
  lower[,i] <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
  upper[,i] <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
  control[,i] <- out_data[[i+20]]$total_M_arabiensis + out_data[[i+20]]$total_M_funestus + out_data[[i+20]]$total_M_gambiae
}
lower <- apply(lower, MARGIN = 1, FUN = mean)
upper <- apply(upper, MARGIN = 1, FUN = mean)
control <- apply(control, MARGIN = 1, FUN = mean)
plot(out_data$Control$timestep/365 +2000, 
     control,
     type="l", frame.plot = F, lwd=2, xlim = c(2019,2022),
     xlab = "Year", ylab = "Total mosquitoes", ylim=c(0,3000000))
polygon(c(out_data$Control$timestep/365 +2000, rev(out_data$Control$timestep/365 +2000)),
        c(lower,rev(upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
        border = F)

plot(out_data$Control$timestep/365 +2000, 
     out_data$Control$n_detect_730_3649/out_data$Control$n_730_3649,
     type="l", frame.plot = F, xlim = c(2019,2022), ylim=c(0,0.8),
     xlab = "Year", ylab = "PfPr2-10")
# i=2
# total_M_lower <- out_data[[i]]$total_M_arabiensis + out_data[[i]]$total_M_funestus + out_data[[i]]$total_M_gambiae
# total_M_upper <- out_data[[i+10]]$total_M_arabiensis + out_data[[i+10]]$total_M_funestus + out_data[[i+10]]$total_M_gambiae
# polygon(c(out_data[[i]]$timestep/365 +2000, rev(out_data[[i+10]]$timestep/365 +2000)),
#         c(total_M_lower,rev(total_M_upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
#         border = F)
for (i in 1:10) {
  pfpr_lower <- out_data[[i]]$n_detect_730_3649/out_data[[i]]$n_730_3649
  pfpr_upper <- out_data[[i+10]]$n_detect_730_3649/out_data[[i+10]]$n_730_3649
  pfpr_control <- out_data[[i+20]]$n_detect_730_3649/out_data[[i+20]]$n_730_3649
  lines(out_data[[i]]$timestep/365+2000, pfpr_lower, col = brewer.pal(10, "RdYlBu")[2])
  lines(out_data[[i+10]]$timestep/365+2000, pfpr_upper, col = brewer.pal(10, "RdYlBu")[9])
  lines(out_data[[i+20]]$timestep/365+2000, pfpr_control, col="black")
}
lower <- matrix(0,nrow = 8395, ncol=10)
upper <- matrix(0,nrow = 8395, ncol=10)
control <- matrix(0,nrow = 8395, ncol=10)
for (i in 1:10) {
  lower[,i] <- out_data[[i]]$n_detect_730_3649/out_data[[i]]$n_730_3649
  upper[,i] <- out_data[[i+10]]$n_detect_730_3649/out_data[[i+10]]$n_730_3649
  control[,i] <- out_data[[i+20]]$n_detect_730_3649/out_data[[i+20]]$n_730_3649
}
lower <- apply(lower, MARGIN = 1, FUN = mean)
upper <- apply(upper, MARGIN = 1, FUN = mean)
control <- apply(control, MARGIN = 1, FUN = mean)
plot(out_data$Control$timestep/365 +2000, 
     control,
     type="l", frame.plot = F, lwd=2, xlim = c(2019,2022),
     xlab = "Year", ylab = "PfPr2-10", ylim=c(0,0.2))
polygon(c(out_data$Control$timestep/365 +2000, rev(out_data$Control$timestep/365 +2000)),
        c(lower,rev(upper)), col = adjustcolor(col = "orange", alpha.f = 0.5),
        border = F)


int <- western_rural$interventions

