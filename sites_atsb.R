################################################################################
## The purpose of this script is to attempt to generate model runs from       ##
## malariasimulation across the three countries involved in the ongoing ATSB  ##
## trials (Mali, Zambia and Kenya). The ideal outcome would be to save output ##
## data from all admin 1 units from 2000 to 2028 against which different      ##
## intervention packages, including ATSBs, could be compared.                 ##
################################################################################
setwd("Documents/GitHub/")
library(foresite)
library(site)
library(devtools)
# install_github("nmoghaddas19/malariasimulation")
library(malariasimulation)
out_mali <- readRDS("OneDrive_1_2-13-2023/out_mali.RDS")
out_zambia <- readRDS("OneDrive_1_2-13-2023/out_zambia.RDS")
out_kenya <- readRDS("OneDrive_1_2-13-2023/out_kenya.RDS")
out_atsb_mali <- readRDS("OneDrive_1_2-13-2023/out_atsb_mali.RDS")

################################################################################
###############################   Site 1: Mali   ###############################
################################################################################

mali <- MLI
names(mali)
mali$eir
mali_sites <- list()
mali_params <- list()
for (i in 1:nrow(mali$sites)) {
  name <- paste0(mali$sites$name_1[i],"_", mali$sites$urban_rural[i])
  mali_sites[[name]] <- single_site(site_file = mali, i)
  params <- site_parameters(
      interventions = mali_sites[[i]]$interventions,
      demography = mali_sites[[i]]$demography,
      vectors = mali_sites[[i]]$vectors,
      seasonality = mali_sites[[i]]$seasonality,
      eir = mali_sites[[i]]$eir$eir[1],
      overrides = list(human_population = 5000)
    )
  params_name <- paste0(name, "_", "params")
  mali_params[[params_name]] <- params
}

par(mfrow=c(3,3), las=1)
out_mali <- list()
for (i in 1:length(mali_params)) {
  name <- paste0(mali$sites$name_1[i],"_", mali$sites$urban_rural[i], "_", "data")
  out_mali[[name]] <- run_simulation(timesteps = mali_params[[i]]$timesteps + 5*365, 
                                parameters = mali_params[[i]])
  plot(out_mali[[i]]$timestep/365 +2000, out_mali[[i]]$n_detect_730_3649/out_mali[[i]]$n_730_3649,
       ylim=c(0,1), frame.plot=F, xlab = "Year", ylab = "PfPr2-10", type = "l", 
       lwd = 2.5, xlim = c(2015, 2030))
  title(names(mali_sites)[i])
}

saveRDS(object = out_mali, file = "out_mali.RDS")

# inclusing atsbs in 2017
mali_atsb_sites <- list()
mali_atsb_params <- list()
for (i in 1:nrow(mali$sites)) {
  name <- paste0(mali$sites$name_1[i],"_", mali$sites$urban_rural[i])
  mali_atsb_sites[[name]] <- single_site(site_file = mali, i)
  params <- site_parameters(
    interventions = mali_sites[[i]]$interventions,
    demography = mali_sites[[i]]$demography,
    vectors = mali_sites[[i]]$vectors,
    seasonality = mali_sites[[i]]$seasonality,
    eir = mali_sites[[i]]$eir$eir[1],
    overrides = list(human_population = 5000,
                     mu_atsb = c(0.35,0.35,0.35))
  )
  params <- set_atsb(parameters = params, timesteps = (17*365):(18*365), coverages = rep(1,366))
  params_name <- paste0(name, "_", "params")
  mali_params[[params_name]] <- params
}

par(mfrow=c(3,3), las=1)
out_atsb_mali <- list()
for (i in 1:length(mali_params)) {
  name <- paste0(mali$sites$name_1[i],"_", mali$sites$urban_rural[i], "_", "data")
  out_atsb_mali[[name]] <- run_simulation(timesteps = mali_params[[i]]$timesteps, 
                                     parameters = mali_params[[i]])
  plot(out_atsb_mali[[i]]$timestep/365 +2000, out_atsb_mali[[i]]$n_detect_730_3649/out_atsb_mali[[i]]$n_730_3649,
       ylim=c(0,1), frame.plot=F, xlab = "Year", ylab = "PfPr2-10", type = "l", col="dodgerblue", 
       lwd = 2.5, xlim = c(2015, 2030))
  title(names(mali_sites)[i])
}

saveRDS(out_atsb_mali, "out_atsb_mali.RDS")

################################################################################
##############################   Site 2: Zambia   ##############################
################################################################################

zambia <- ZMB
zambia_sites <- list()
zambia_params <- list()
for (i in 1:nrow(zambia$sites)) {
  name <- paste0(zambia$sites$name_1[i],"_", zambia$sites$urban_rural[i])
  zambia_sites[[name]] <- single_site(site_file = zambia, i)
  params <- site_parameters(
    interventions = zambia_sites[[i]]$interventions,
    demography = zambia_sites[[i]]$demography,
    vectors = zambia_sites[[i]]$vectors,
    seasonality = zambia_sites[[i]]$seasonality,
    eir = zambia_sites[[i]]$eir$eir[1],
    overrides = list(human_population = 5000
                     # mu_atsb = c(0.039, 0.089, 0.039)
                     )
  )
  #params <- set_atsb(parameters = params, timesteps = (23*365):(24*365), coverages = rep(1,366))
  params_name <- paste0(name, "_", "params")
  zambia_params[[params_name]] <- params
}

plot.new()
par(mfrow=c(3,3), las=1)
out_zambia <- list()
for (i in 1:length(zambia_params)) {
  name <- paste0(zambia$sites$name_1[i],"_", zambia$sites$urban_rural[i], "_", "data")
  out_zambia[[name]] <- run_simulation(timesteps = zambia_params[[i]]$timesteps + 5*365, 
                                 parameters = zambia_params[[i]])
  plot(out_zambia[[i]]$timestep/365 +2000, out_zambia[[i]]$n_detect_730_3649/out_zambia[[i]]$n_730_3649,
       ylim=c(0,1), frame.plot=F, xlab = "Year", ylab = "PfPr2-10", type = "l", 
       lwd = 2.5, xlim = c(2015, 2030))
  title(names(zambia_sites)[i])
}

saveRDS(out_zambia, "out_zambia.RDS")

################################################################################
##############################   Site 3: Kenya   ###############################
################################################################################

kenya <- KEN
kenya_sites <- list()
kenya_params <- list()
for (i in 1:nrow(kenya$sites)) {
  name <- paste0(kenya$sites$name_1[i],"_", kenya$sites$urban_rural[i])
  kenya_sites[[name]] <- single_site(site_file = kenya, i)
  params <- site_parameters(
    interventions = kenya_sites[[i]]$interventions,
    demography = kenya_sites[[i]]$demography,
    vectors = kenya_sites[[i]]$vectors,
    seasonality = kenya_sites[[i]]$seasonality,
    eir = kenya_sites[[i]]$eir$eir[1],
    overrides = list(human_population = 5000)
  )
  params_name <- paste0(name, "_", "params")
  kenya_params[[params_name]] <- params
}

plot.new()
par(mfrow=c(3,3), las=1)
out_kenya <- list()
for (i in 1:length(kenya_params)) {
  name <- paste0(kenya$sites$name_1[i],"_", kenya$sites$urban_rural[i], "_", "data")
  out_kenya[[name]] <- run_simulation(timesteps = kenya_params[[i]]$timesteps, 
                                parameters = kenya_params[[i]])
  plot(out[[i]]$timestep/365 +2000, out[[i]]$n_detect_730_3649/out[[i]]$n_730_3649,
       ylim=c(0,1), frame.plot=F, xlab = "Year", ylab = "PfPr2-10", type = "l", 
       lwd = 2.5, xlim = c(2000, 2025))
  title(names(kenya_sites)[i])
}

saveRDS(out_kenya, "out_kenya.RDS")
  
  
################################################################################
## The aim here is to use the ATSB implemented version of malariasimulation   ##
## that I have developed to generate meaningful confidence intervals from the ##
## available data. My approach here will be twofold. Firstly, I will use the  ##
## uncertainty in the count data to find what the corresponding uncertainty   ##
## in the bait feeding rate would be. Secondly, I will use the uncertainty in ##
## the feeding rate data to estimate what uncertainty in the count data we    ##
## would expect. Therefore, I will write a script that takes each of those    ##
## uncertainties and generates the other. I will use the site files for rural ## 
## western Mali which is where the trial took place.                          ##
################################################################################

plot(out_mali$Bamako_rural_data$timestep/365 + 2000, 
     out_mali$Bamako_rural_data$n_detect_730_3649/out_mali$Kayes_rural_data$n_730_3649,
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", ylab = "PfPr2-10", 
     xlim = c(2010,2020), ylim = c(0,1))
# title("")
plot(out_mali$Kayes_rural_data$timestep/365 + 2000, 
     out_mali$Kayes_rural_data$n_detect_730_3649/out_mali$Kayes_rural_data$n_730_3649,
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", ylab = "PfPr2-10", 
     xlim = c(2010,2020), ylim = c(0,1))
lines(out_atsb_mali$Kayes_rural_data$timestep/365 + 2000,
      out_atsb_mali$Kayes_rural_data$n_detect_730_3649/out_atsb_mali$Kayes_rural_data$n_730_3649,
      lwd = 2, col = "dodgerblue")
plot(out_mali$Kayes_rural_data$timestep/365 + 2000, 
     out_mali$Kayes_rural_data[["Im_gambiae_count"]] + 
       out_mali$Kayes_rural_data[["Sm_gambiae_count"]] + 
       out_mali$Kayes_rural_data[["Pm_gambiae_count"]],
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", ylab = "An. gambiae count", 
     xlim = c(2010,2020))
lines(out_atsb_mali$Kayes_rural_data$timestep/365 + 2000, 
      out_atsb_mali$Kayes_rural_data[["Im_gambiae_count"]] + 
        out_atsb_mali$Kayes_rural_data[["Sm_gambiae_count"]] + 
        out_atsb_mali$Kayes_rural_data[["Pm_gambiae_count"]], 
      lwd = 2, col = "dodgerblue")

feeding_rate <- c(0.35, 0.35, 0.35) # arab, fun, gamb
mali <- MLI
mali$sites
kayes_rural <- single_site(mali, 5)
kayes_rural_params_lower <- site_parameters(
  interventions = kayes_rural$interventions,
  demography = kayes_rural$demography,
  vectors = kayes_rural$vectors,
  seasonality = kayes_rural$seasonality,
  eir = kayes_rural$eir$eir[1],
  overrides = list(human_population = 5000,
                   mu_atsb = c(0.15,0.15,0.15))
)
kayes_rural_params_lower <- set_atsb(parameters = kayes_rural_params_lower,
                               timesteps = (17*365+5*30):(18*365), 
                               coverages = rep(1,366-5*30))
out_kayes_rural_lower <- run_simulation(timesteps = kayes_rural_params_lower$timesteps,
                                  parameters = kayes_rural_params_lower)

kayes_rural_params_upper <- site_parameters(
  interventions = kayes_rural$interventions,
  demography = kayes_rural$demography,
  vectors = kayes_rural$vectors,
  seasonality = kayes_rural$seasonality,
  eir = kayes_rural$eir$eir[1],
  overrides = list(human_population = 5000,
                   mu_atsb = c(0.17,0.17,0.17))
)
kayes_rural_params_upper <- set_atsb(parameters = kayes_rural_params_upper,
                                     timesteps = (17*365+5*30):(18*365), 
                                     coverages = rep(1,366-5*30))
out_kayes_rural_upper <- run_simulation(timesteps = kayes_rural_params_upper$timesteps,
                                        parameters = kayes_rural_params_upper)

plot(out_mali$Kayes_rural_data$timestep/365 + 2000, 
     out_mali$Kayes_rural_data$n_detect_730_3649/out_mali$Kayes_rural_data$n_730_3649,
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", ylab = "PfPr2-10", 
     xlim = c(2016,2018), ylim = c(0,1))
lines(out_kayes_rural_upper$timestep/365 + 2000,
      out_kayes_rural_upper$n_detect_730_3649/out_kayes_rural_upper$n_730_3649,
      lwd = 2, col = "dodgerblue")
lines(out_kayes_rural_lower$timestep/365 + 2000,
      out_kayes_rural_lower$n_detect_730_3649/out_kayes_rural_upper$n_730_3649,
      lwd = 2, col = "royalblue")
polygon(c(out_kayes_rural_upper$timestep/365 + 2000, rev(out_kayes_rural_lower$timestep/365 + 2000)),
        c(out_kayes_rural_upper$n_detect_730_3649/out_kayes_rural_upper$n_730_3649, 
          rev(out_kayes_rural_lower$n_detect_730_3649/out_kayes_rural_lower$n_730_3649)),
        col=adjustcolor("royalblue", alpha.f=0.30), border=FALSE)
        
feed <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)
o <- list()
for (i in 1:10) {
  name <- as.character(feed[i])
  p <- site_parameters(
    interventions = kayes_rural$interventions,
    demography = kayes_rural$demography,
    vectors = kayes_rural$vectors,
    seasonality = kayes_rural$seasonality,
    eir = kayes_rural$eir$eir[1],
    overrides = list(human_population = 5000,
                     mu_atsb = c(feed[i],feed[i],feed[i]))
  )
  p <- set_atsb(parameters = p,
           timesteps = (17*365+5*30):(18*365), 
           coverages = rep(1,366-5*30))
  o[[name]] <- run_simulation(timesteps = p$timesteps,
                              parameters = p)
  print(i)
}
display.brewer.all(type="div")
brewer.pal(10, "RdYlBu")
plot(out_kayes_rural$timestep/365 + 2000, 
     out_kayes_rural$Sm_gambiae_count +
       out_kayes_rural$Pm_gambiae_count +
       out_kayes_rural$Im_gambiae_count,
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", 
     xlim = c(2016,2018), ylab = "An. gambiae population")
for (i in 1:9) {
  polygon(c(o[[i]]$timestep/365 + 2000, rev(o[[i+1]]$timestep/365 + 2000)),
          c(o[[i]]$Sm_gambiae_count+o[[i]]$Pm_gambiae_count+o[[i]]$Im_gambiae_count,
            rev(o[[i+1]]$Sm_gambiae_count+o[[i+1]]$Pm_gambiae_count+o[[i+1]]$Im_gambiae_count)), 
          col=adjustcolor(brewer.pal(9, "RdYlBu")[i], alpha.f=0.50), border=FALSE)
}
legend(x="topright", legend = c(as.character(feed)[2:10]), 
       fill = adjustcolor(brewer.pal(9, "RdYlBu"), alpha.f=0.50),
       bty = "n")

plot(out_mali$Kayes_rural_data$timestep/365 + 2000, 
     out_mali$Kayes_rural_data$n_detect_730_3649/out_mali$Kayes_rural_data$n_730_3649,
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", ylab = "PfPr2-10", 
     xlim = c(2016,2018), ylim = c(0,1))
for (i in 1:9) {
  polygon(c(o[[i]]$timestep/365 + 2000, rev(o[[i+1]]$timestep/365 + 2000)),
          c(o[[i]]$n_detect_730_3649/o[[i]]$n_730_3649,
            rev(o[[i+1]]$n_detect_730_3649/o[[i+1]]$n_730_3649)), 
          col=adjustcolor(brewer.pal(9, "RdYlBu")[i], alpha.f=0.50), border=FALSE)
}
legend(x="topright", legend = c(as.character(feed)[2:10]), 
       fill = adjustcolor(brewer.pal(9, "RdYlBu"), alpha.f=0.50),
       bty = "n")
plot(o[[1]]$timestep/365 + 2000, (o[[1]]$total_M_gambiae-o[[2]]$total_M_gambiae)/o[[1]]$total_M_gambiae, 
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", ylab = "% reduction", 
     xlim = c(2016,2018), ylim = c(0,1))
peak_season_offset(parameters= p)

# the purpose of this code is to visualise the confidence intervals quoted by 
# Fraser et al. The quoted mean reduction in mosquito catch was 57% [33-72% 95%CI]
peak_season_offset(parameters= p)
# it seems like 40% feeding rate gives a 73% reduction - fairly close to 72
(o[[1]]$total_M_gambiae[17*365+228]-o[[4]]$total_M_gambiae[17*365+228])/
  o[[1]]$total_M_gambiae[17*365+228]
# 10% feeding rate gives 38% reduction but 5% feeding gives 23% so probably need 8 or 9%
(o[[1]]$total_M_gambiae[17*365+228]-o[[10]]$total_M_gambiae[17*365+228])/
  o[[1]]$total_M_gambiae[17*365+228]
#lets try 8% feeding rate
p <- site_parameters(
  interventions = kayes_rural$interventions,
  demography = kayes_rural$demography,
  vectors = kayes_rural$vectors,
  seasonality = kayes_rural$seasonality,
  eir = kayes_rural$eir$eir[1],
  overrides = list(human_population = 5000,
                   mu_atsb = c(0.105,0.105,0.105))
)
p <- set_atsb(parameters = p,
              timesteps = (17*365+5*30):(18*365), 
              coverages = rep(1,366-5*30))
o10.5 <- run_simulation(timesteps = p$timesteps,
               parameters = p)
# 7.8% gives the lower bound of 33%
(o[[1]]$total_M_gambiae[17*365+228]-o10.5$total_M_gambiae[17*365+228])/
  o[[1]]$total_M_gambiae[17*365+228]

p <- site_parameters(
  interventions = kayes_rural$interventions,
  demography = kayes_rural$demography,
  vectors = kayes_rural$vectors,
  seasonality = kayes_rural$seasonality,
  eir = kayes_rural$eir$eir[1],
  overrides = list(human_population = 5000,
                   mu_atsb = c(0.38,0.38,0.38))
)
p <- set_atsb(parameters = p,
              timesteps = (17*365+5*30):(18*365), 
              coverages = rep(1,366-5*30))
o38 <- run_simulation(timesteps = p$timesteps,
                       parameters = p)
# 38% gives the upper bound of 72%
(o[[1]]$total_M_gambiae[17*365+228]-o50$total_M_gambiae[17*365+228])/
  o[[1]]$total_M_gambiae[17*365+228]

plot(out_kayes_rural$timestep/365 + 2000, 
     out_kayes_rural$Sm_gambiae_count +
       out_kayes_rural$Pm_gambiae_count +
       out_kayes_rural$Im_gambiae_count,
     type = "l", lwd = 2, frame.plot = F, xlab = "Year", 
     xlim = c(2016,2018), ylab = "An. gambiae population")
polygon(c(o7.8$timestep/365 + 2000, rev(o38$timestep/365 +2000)),
        c(o7.8$total_M_gambiae, rev(o38$total_M_gambiae)),
        col= adjustcolor("dodgerblue", alpha.f=0.50), border = F)
# what does that sugar feeding uncertainty look like on the measured data
sugar_feeding <- read.csv("DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
sugar_feeding <- sugar_feeding[-(64:65),]
par(las=1)
plot.default(sugar_feeding$month, 
             sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2,
             cex=0.9, pch=20, frame.plot = F, xlab = "Month", ylab = "% bait fed",
             ylim = c(0,1))
qt(c(0.025,0.975), df = 6)
sugar_feeding |>
  group_by(month) |>
  summarise(mean = mean(females.ASB.positive/TOTAL.Sample.females.Day.2),
            CI = 2.44*sd(females.ASB.positive/TOTAL.Sample.females.Day.2)/sqrt(n())) -> sugar_feeding_grouped
polygon(c(4,12,12,4), c(0.15,0.15,0.45,0.45),
        col = adjustcolor("grey", alpha.f = 0.5), border = F)
lines(sugar_feeding_grouped$month, sugar_feeding_grouped$mean, col="red", lwd=2)
arrows(sugar_feeding_grouped$month, 
       sugar_feeding_grouped$mean-sugar_feeding_grouped$CI,
       sugar_feeding_grouped$month,
       sugar_feeding_grouped$mean+sugar_feeding_grouped$CI,
       angle = 90,
       code = 3,
       length = 0.1, 
       col="red")
# okay what does the uncertainty in feeding rate look like overlayed on the 
# count data

# step 1: load in the count data for both years
hlc_may_jul <- read.csv("HLC data 2016 (unchanged) processed/HLC May_Jul2016-Table 1.csv")[,-13]
hlc_aug_sep <- read.csv("HLC data 2016 (unchanged) processed/HLC Aug_Sept 2016-Table 1.csv")[,-13]
hlc_oct_dec <- read.csv("HLC data 2016 (unchanged) processed/HLC Oct_Dber-Table 1.csv")[,-13]
hlc_aug_sep <- hlc_aug_sep[,-16]
hlc_apr_jul_2017 <- read.csv("HLC data processed/April_May_June_July-Table 1.csv")[,1:16][,-13]
hlc_aug_sep_2017 <- read.csv("HLC data processed/August_September-Table 1.csv")[,1:16][,-13]
hlc_oct_dec_2017 <- read.csv("HLC data processed/Oct_Nvber_Dcber-Table 1.csv")[,1:16][,-13]
hlc <- rbind(hlc_may_jul, hlc_aug_sep, hlc_oct_dec,
             hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)
rm(hlc_may_jul, hlc_aug_sep, hlc_oct_dec,
   hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)

table(hlc$Study.Sites)
hlc <- hlc[-which(hlc$Position=="N/A"),]

hlc |>
  mutate(time_frame = ifelse(NightFrames=="18h-20h", 1,
                             ifelse(NightFrames=="20-22h"|NightFrames=="20h-22h", 2,
                                    ifelse(NightFrames=="22h-00h"|NightFrames=="22H-00h", 3,
                                           ifelse(NightFrames=="00h-02h", 4,
                                                  ifelse(NightFrames=="02h-04h"|NightFrames=="02h-4h", 5,
                                                         ifelse(NightFrames=="04h-06h"|NightFrames=="04H-06h", 6,
                                                                NightFrames))))))) |> 
  mutate(Site = ifelse(Study.Sites=="Kignélé"|Study.Sites=="Kignele ", "Kignele",
                       ifelse(Study.Sites=="Niaganabougou"|Study.Sites=="Nianganabougou", "Nianguanabougou",
                              ifelse(Study.Sites=="Sirakélé"|Study.Sites=="sirakele","Sirakele",
                                     ifelse(Study.Sites=="Krékrélo","Krekrelo",
                                            ifelse(Study.Sites=="Trekourou", "Trekrou",
                                                   Study.Sites)))))) |>
  mutate(Treatment = ifelse(Site %in% c("Krekrelo", "Sirakele", "Trekrou", "Farabale", "Kignele", "Tiko", "Sambadani"), 1, 0)) |>
  mutate(Position = ifelse(Position=="Indoor", "Inside",
                           ifelse(Position=="Outdoor"|Position=="outside", "Outside",
                                  Position))) |>
  mutate(Month = ifelse(Month=="Août", "August",
                        ifelse(Month=="Juillet"|Month=="Jul","July",
                               ifelse(Month=="Jun", "June", 
                                      ifelse(Month=="S+3835:3869eptember", "September",
                                             Month))))) |>
  group_by(Year, Month, Site) |>
  summarise(Count=n(), Outdoor=sum(Position=="Outside"), Indoor=sum(Position=="Inside")) -> hlc_bysitemonth
hlc_bysitemonth$Month_Year <- paste0(hlc_bysitemonth$Month,hlc_bysitemonth$Year)
hlc_bysitemonth |>
  mutate(Treatment = ifelse(Site %in% c("Krekrelo", "Sirakele", "Trekrou", "Farabale", "Kignele", "Tiko", "Sambadani"), 1, 0)) |>
  mutate(Month_Year = parse_date_time(Month_Year, "my")) -> hlc_bysitemonth
hlc_bysitemonth |>
  group_by()
plot(hlc_bysitemonth$Month_Year, hlc_bysitemonth$Count, pch=20, frame.plot = F,
     col=hlc_bysitemonth$Treatment+1, ylab = "Count", xlab=NA)
#
glmer(
  cbind(females.ASB.positive, TOTAL.Sample.females.Day.2 - females.ASB.positive) ~ month + (1|Village),
  family = binomial, data = sugar_feeding) |>
  summary()
InvLogit(-0.1694-0.2221+0.3349*c(1.96, 0, -1.96))



        