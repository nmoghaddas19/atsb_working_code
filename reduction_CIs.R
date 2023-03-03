# use bootstrap to estimate percent reduction in mean count and uncertainty around it
# cdc
# import data
cdc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/CDC-Table 1.csv")

# bootstrap
mean_reduction <- c()
quantiles_reduction <- matrix(0, nrow=9, ncol=3)
for (j in 4:12) {
  cdc_2017 |>
    filter(Month==j & Experimental.or.control=="Con.") -> control
  cdc_2017 |>
    filter(Month==j & Experimental.or.control=="Exp.") -> experiment
  for (i in 1:5000) {
    control_sample <- sample(control$tot.f, size = 7, replace = T)
    experiment_sample <- sample(experiment$tot.f, size = 7, replace = T)
    mean_reduction[i] <- (mean(control_sample)-mean(experiment_sample)) / mean(control_sample)
  }
  quantiles_reduction[j-3,] <- quantile(mean_reduction, prob=c(0.025,0.5,0.975))
}
# get right days_shifted from shift.R
plot(2017+4:12/12, 1-quantiles_reduction[,2], type="l", lwd=2, frame.plot = F, ylim=c(0,1),
     xlab="Year", ylab="Proportion of control population", xlim=c(2017.6,2018))
polygon(c(2017+4:12/12, rev(2017+4:12/12)),
        c(1-quantiles_reduction[,1], rev(1-quantiles_reduction[,3])),
        col = adjustcolor("black", alpha.f = 0.1), border = FALSE)
q <- out_bounds[[1]]$total_M_gambiae/malariasim_control$total_M_gambiae
w <- out_bounds[[2]]$total_M_gambiae/malariasim_control$total_M_gambiae
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(q, 
          rev(w)),
        col = adjustcolor("dodgerblue", alpha.f = 0.4), border = FALSE)
title("CDC traps")

# malaise
# import data
malaise_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/Malaise-Table 1.csv")

# bootstrap
mean_reduction <- c()
quantiles_reduction <- matrix(0, nrow=9, ncol=3)
for (j in 4:12) {
  malaise_2017 |>
    filter(Month==j & Experimental.or.control=="Con.") -> control
  malaise_2017 |>
    filter(Month==j & Experimental.or.control=="Exp.") -> experiment
  for (i in 1:5000) {
    control_sample <- sample(control$tot.f, size = 7, replace = T)
    experiment_sample <- sample(experiment$tot.f, size = 7, replace = T)
    mean_reduction[i] <- (mean(control_sample)-mean(experiment_sample)) / mean(control_sample)
  }
  quantiles_reduction[j-3,] <- quantile(mean_reduction, prob=c(0.025,0.5,0.975))
}
# get right days_shifted from shift.R
plot(2017+4:12/12, 1-quantiles_reduction[,2], type="l", lwd=2, frame.plot = F, ylim=c(0,1),
     xlab="Year", ylab="Proportion of control population", xlim=c(2017.6,2018))
polygon(c(2017+4:12/12, rev(2017+4:12/12)),
        c(1-quantiles_reduction[,1], rev(1-quantiles_reduction[,3])),
        col = adjustcolor("black", alpha.f = 0.1), border = FALSE)
q <- out_bounds[[1]]$total_M_gambiae/malariasim_control$total_M_gambiae
w <- out_bounds[[2]]$total_M_gambiae/malariasim_control$total_M_gambiae
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(q, 
          rev(w)),
        col = adjustcolor("dodgerblue", alpha.f = 0.4), border = FALSE)
title("Malaise traps")

# pyrethrum spray catch
# import data
psc_2017 <- read.csv("atsb_working_code/DB CDC Malaise PSC catches 2017/PSC-Table 1.csv")
psc_2017 |>
  filter(!is.na(Month)) -> psc_2017
# bootstrap
mean_reduction <- c()
quantiles_reduction <- matrix(0, nrow=9, ncol=3)
for (j in 4:12) {
  psc_2017 |>
    filter(Month==j & Experimental.or.control=="Con.") -> control
  psc_2017 |>
    filter(Month==j & Experimental.or.control=="Exp.") -> experiment
  for (i in 1:5000) {
    control_sample <- sample(control$tot.f, size = 7, replace = T)
    experiment_sample <- sample(experiment$tot.f, size = 7, replace = T)
    mean_reduction[i] <- (mean(control_sample)-mean(experiment_sample)) / mean(control_sample)
  }
  quantiles_reduction[j-3,] <- quantile(mean_reduction, prob=c(0.025,0.5,0.975))
}
# get right days_shifted from shift.R
plot(2017+4:12/12, 1-quantiles_reduction[,2], type="l", lwd=2, frame.plot = F, ylim=c(0,1),
     xlab="Year", ylab="Proportion of control population", xlim=c(2017.6,2018))
polygon(c(2017+4:12/12, rev(2017+4:12/12)),
        c(1-quantiles_reduction[,1], rev(1-quantiles_reduction[,3])),
        col = adjustcolor("black", alpha.f = 0.1), border = FALSE)
q <- out_bounds[[1]]$total_M_gambiae/malariasim_control$total_M_gambiae
w <- out_bounds[[2]]$total_M_gambiae/malariasim_control$total_M_gambiae
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(q, 
          rev(w)),
        col = adjustcolor("dodgerblue", alpha.f = 0.4), border = FALSE)
title("Pyrethrum spray catch")

# human landing catch
# import data
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
hlc_2017$Month_Year <- 
hlc_bysitemonth$Month <- month(hlc_bysitemonth$Month_Year)
hlc_2017 <- hlc_bysitemonth
# ----
# bootstrap
mean_reduction <- c()
quantiles_reduction <- matrix(0, nrow=9, ncol=3)
for (j in 4:12) {
  hlc_2017 |>
    filter(Month==j & Treatment==0) -> control
  hlc_2017 |>
    filter(Month==j & Treatment==1) -> experiment
  for (i in 1:5000) {
    control_sample <- sample(control$tot.f, size = 7, replace = T)
    experiment_sample <- sample(experiment$tot.f, size = 7, replace = T)
    mean_reduction[i] <- (mean(control_sample)-mean(experiment_sample)) / mean(control_sample)
  }
  quantiles_reduction[j-3,] <- quantile(mean_reduction, prob=c(0.025,0.5,0.975))
}
# get right days_shifted from shift.R
plot(2017+4:12/12, 1-quantiles_reduction[,2], type="l", lwd=2, frame.plot = F, ylim=c(0,1),
     xlab="Year", ylab="Proportion of control population", xlim=c(2017.6,2018))
polygon(c(2017+4:12/12, rev(2017+4:12/12)),
        c(1-quantiles_reduction[,1], rev(1-quantiles_reduction[,3])),
        col = adjustcolor("black", alpha.f = 0.1), border = FALSE)
q <- out_bounds[[1]]$total_M_gambiae/malariasim_control$total_M_gambiae
w <- out_bounds[[2]]$total_M_gambiae/malariasim_control$total_M_gambiae
polygon(c(out_bounds[[1]]$timestep/365+2000+days_shifted/365, rev(out_bounds[[2]]$timestep/365+2000+days_shifted/365)),
        c(q, 
          rev(w)),
        col = adjustcolor("dodgerblue", alpha.f = 0.4), border = FALSE)
title("Pyrethrum spray catch")

