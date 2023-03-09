# The purpose of this script is to write program that will take malariasim
# time series data and shift it to best fit a set of true data

# import standard malariasim run 
malariasim_control <- readRDS("OneDrive_1_2-13-2023/out_mali.RDS")$Kayes_rural_data

# import data
cdc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/CDC-Table 1.csv")


# find scaling factor 
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> cdc_2017_con
scaler <- max(cdc_2017_con$Mean)/max(malariasim_control$total_M_gambiae[(17*365):(18*365)])

# visualise
plot(malariasim_control$timestep/365+2000, malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018))
lines(2017+(cdc_2017_con$Month)/12, cdc_2017_con$Mean, lwd=2, col=1, lty=2)

# calculate sum of least squares
sum_of_squares <- c()

for (i in 1:80) {
  true_values <- cdc_2017_con$Mean 
  y <- malariasim_control$total_M_gambiae*scaler
  x <- malariasim_control$timestep/365 + 2000 + i/365
  model_values <- approx(x, y, xout= 2017 + c(4:12)/12)
  sum_of_squares[i] <- sum((true_values - model_values$y)^2)
}
plot(1:80, sum_of_squares, type="l", lwd=2, frame.plot = F, xlab="Days shifted")

days_shifted <- which(sum_of_squares == min(sum_of_squares))

# lets see the data with that shift
plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018),
     xlab="Year", ylab="Population")
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(out_bounds[[1]]$total_M_gambiae*scaler, rev(out_bounds[[2]]$total_M_gambiae)*scaler),
        col = adjustcolor("dodgerblue", alpha.f = 0.5), border = FALSE)
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> cdc_2017_con
lines(2017+(cdc_2017_con$Month)/12, cdc_2017_con$Mean, lwd=2, col=1, lty=2)
arrows(2017+(cdc_2017_con$Month)/12,
       cdc_2017_con$Mean-cdc_2017_con$CI,
       2017+(cdc_2017_con$Month)/12,
       cdc_2017_con$Mean+cdc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.05,
       col=1)
cdc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> cdc_2017_exp
lines(2017+(cdc_2017_exp$Month)/12, cdc_2017_exp$Mean, lwd=2, col=2, lty=2)
arrows(2017+(cdc_2017_exp$Month)/12, 
       cdc_2017_exp$Mean-cdc_2017_exp$CI,
       2017+(cdc_2017_exp$Month)/12,
       cdc_2017_exp$Mean+cdc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=2)
legend(x="topleft", legend=c("Control data", "ATSB data", "Control model", "ATSB model"), 
       col=c(1,2,1,"dodgerblue"), lwd=2, lty=c(2,2,1,1), bty="n")
title("CDC traps")

# rinse and repeat for the other traps 
# malaise
malaise_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/Malaise-Table 1.csv")

malaise_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> malaise_2017_con
scaler <- max(malaise_2017_con$Mean)/max(malariasim_control$total_M_gambiae[(17*365):(18*365)])

plot(malariasim_control$timestep/365+2000, malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,800), xlim=c(2016,2018))
lines(2017+(malaise_2017_con$Month)/12, malaise_2017_con$Mean, lwd=2, col=1, lty=2)

sum_of_squares <- c()

for (i in 1:80) {
  true_values <- malaise_2017_con$Mean 
  y <- malariasim_control$total_M_gambiae*scaler
  x <- malariasim_control$timestep/365 + 2000 + i/365
  model_values <- approx(x, y, xout= 2017 + c(4:12)/12)
  sum_of_squares[i] <- sum((true_values - model_values$y)^2)
}
plot(1:80, sum_of_squares, type="l", lwd=2, frame.plot = F, xlab="Days shifted")

days_shifted <- which(sum_of_squares == min(sum_of_squares))

plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,600), xlim=c(2016,2018),
     xlab="Year", ylab="Population")
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(out_bounds[[1]]$total_M_gambiae*scaler, rev(out_bounds[[2]]$total_M_gambiae)*scaler),
        col = adjustcolor("dodgerblue", alpha.f = 0.5), border = FALSE)
malaise_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> malaise_2017_con
lines(2017+malaise_2017_con$Month/12, malaise_2017_con$Mean, lwd=2, lty=2, col=1)
arrows(2017+malaise_2017_con$Month/12, 
       malaise_2017_con$Mean-malaise_2017_con$CI,
       2017+malaise_2017_con$Month/12,
       malaise_2017_con$Mean+malaise_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=1)
malaise_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> malaise_2017_exp
lines(2017+malaise_2017_exp$Month/12, malaise_2017_exp$Mean, lwd=2, lty=2, col=2)
arrows(2017+malaise_2017_exp$Month/12, 
       malaise_2017_exp$Mean-malaise_2017_exp$CI,
       2017+malaise_2017_exp$Month/12,
       malaise_2017_exp$Mean+malaise_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=2)
legend(x="topleft", legend=c("Control data", "ATSB data", "Control model", "ATSB model"), 
       col=c(1,2,1,"dodgerblue"), lwd=2, lty=c(2,2,1,1), bty="n")
title("Malaise traps")

# pyrethrum spray catch
psc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/PSC-Table 1.csv")
psc_2017 |>
  filter(!is.na(Month)) -> psc_2017

psc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> psc_2017_con
scaler <- max(psc_2017_con$Mean)/max(malariasim_control$total_M_gambiae[(17*365):(18*365)])

plot(malariasim_control$timestep/365+2000, malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,400), xlim=c(2016,2018))
lines(2017+(psc_2017_con$Month)/12, psc_2017_con$Mean, lwd=2, col=1, lty=2)

sum_of_squares <- c()

for (i in 1:80) {
  true_values <- psc_2017_con$Mean 
  y <- malariasim_control$total_M_gambiae*scaler
  x <- malariasim_control$timestep/365 + 2000 + i/365
  model_values <- approx(x, y, xout= 2017 + c(4:12)/12)
  sum_of_squares[i] <- sum((true_values - model_values$y)^2)
}
plot(1:80, sum_of_squares, type="l", lwd=2, frame.plot = F, xlab="Days shifted")

days_shifted <- which(sum_of_squares == min(sum_of_squares))

plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,400), xlim=c(2016,2018),
     xlab="Year", ylab="Population")
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(out_bounds[[1]]$total_M_gambiae*scaler, rev(out_bounds[[2]]$total_M_gambiae)*scaler),
        col = adjustcolor("dodgerblue", alpha.f = 0.5), border = FALSE)
psc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Con.") -> psc_2017_con
lines(2017+psc_2017_con$Month/12, psc_2017_con$Mean, lwd=2, lty=2, col=1)
arrows(2017+psc_2017_con$Month/12, 
       psc_2017_con$Mean-psc_2017_con$CI,
       2017+psc_2017_con$Month/12,
       psc_2017_con$Mean+psc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=1)
psc_2017 |>
  group_by(Month, Experimental.or.control) |>
  summarise(Mean=mean(tot.f),CI=1.96*sd(tot.f)/sqrt(n())) |>
  filter(Experimental.or.control=="Exp.") -> psc_2017_exp
lines(2017+psc_2017_exp$Month/12, psc_2017_exp$Mean, lwd=2, lty=2, col=2)
arrows(2017+psc_2017_exp$Month/12, 
       psc_2017_exp$Mean-psc_2017_exp$CI,
       2017+psc_2017_exp$Month/12,
       psc_2017_exp$Mean+psc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=2)
legend(x="topleft", legend=c("Control data", "ATSB data", "Control model", "ATSB model"), 
       col=c(1,2,1,"dodgerblue"), lwd=2, lty=c(2,2,1,1), bty="n")
title("Pyrethrum spray catch")

# human landing catch
hlc_apr_jul_2017 <- read.csv("atsb_working_code/HLC data processed/April_May_June_July-Table 1.csv")[,1:16][,-13]
hlc_aug_sep_2017 <- read.csv("atsb_working_code/HLC data processed/August_September-Table 1.csv")[,1:16][,-13]
hlc_oct_dec_2017 <- read.csv("atsb_working_code/HLC data processed/Oct_Nvber_Dcber-Table 1.csv")[,1:16][,-13]
hlc <- rbind(hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)
rm(hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)
# pre processing ----
table(hlc$Month)
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
hlc_bysitemonth$Month <- month(hlc_bysitemonth$Month_Year)
hlc_2017 <- hlc_bysitemonth


# ----
hlc_2017 |>
  group_by(Month, Treatment) |>
  summarise(Mean=mean(Count),CI=1.96*sd(Count)/sqrt(n())) |>
  filter(Treatment==0) -> hlc_2017_con
scaler <- max(hlc_2017_con$Mean)/max(malariasim_control$total_M_gambiae[(17*365):(18*365)])

plot(malariasim_control$timestep/365+2000, malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,600), xlim=c(2016,2018))
lines(2017+(hlc_2017_con$Month)/12, hlc_2017_con$Mean, lwd=2, col=1, lty=2)

sum_of_squares <- c()

for (i in 1:80) {
  true_values <- hlc_2017_con$Mean 
  y <- malariasim_control$total_M_gambiae*scaler
  x <- malariasim_control$timestep/365 + 2000 + i/365
  model_values <- approx(x, y, xout= 2017 + c(4:12)/12)
  sum_of_squares[i] <- sum((true_values - model_values$y)^2)
}
plot(1:80, sum_of_squares, type="l", lwd=2, frame.plot = F, xlab="Days shifted")

days_shifted <- which(sum_of_squares == min(sum_of_squares))

plot(malariasim_control$timestep/365+2000 + days_shifted/365, 
     malariasim_control$total_M_gambiae*scaler,
     type="l", lwd=2, frame.plot = F, ylim=c(0,600), xlim=c(2016,2018),
     xlab="Year", ylab="Population")
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(out_bounds[[1]]$total_M_gambiae*scaler, rev(out_bounds[[2]]$total_M_gambiae)*scaler),
        col = adjustcolor("dodgerblue", alpha.f = 0.5), border = FALSE)
hlc_2017 |>
  group_by(Month, Treatment) |>
  summarise(Mean=mean(Count),CI=1.96*sd(Count)/sqrt(n())) |>
  filter(Treatment==0) -> hlc_2017_con
lines(2017+hlc_2017_con$Month/12, hlc_2017_con$Mean, lwd=2, lty=2, col=1)
arrows(2017+hlc_2017_con$Month/12, 
       hlc_2017_con$Mean-hlc_2017_con$CI,
       2017+hlc_2017_con$Month/12,
       hlc_2017_con$Mean+hlc_2017_con$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=1)
hlc_2017 |>
  group_by(Month, Treatment) |>
  summarise(Mean=mean(Count),CI=1.96*sd(Count)/sqrt(n())) |>
  filter(Treatment==1) -> hlc_2017_exp
lines(2017+hlc_2017_exp$Month/12, hlc_2017_exp$Mean, lwd=2, lty=2, col=2)
arrows(2017+hlc_2017_exp$Month/12, 
       hlc_2017_exp$Mean-hlc_2017_exp$CI,
       2017+hlc_2017_exp$Month/12,
       hlc_2017_exp$Mean+hlc_2017_exp$CI,
       angle = 90,
       code = 3,
       length = 0.05, 
       col=2)
legend(x="topleft", legend=c("Control data", "ATSB data", "Control model", "ATSB model"), 
       col=c(1,2,1,"dodgerblue"), lwd=2, lty=c(2,2,1,1), bty="n")
title("Human landing catch")

