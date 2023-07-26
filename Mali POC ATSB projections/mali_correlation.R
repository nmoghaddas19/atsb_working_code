mali <- read.csv("~/Documents/GitHub/atsb_working_code/DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")[-(64:65),1:16]
mali$dyed_fraction <- mali$females.ASB.positive/mali$TOTAL.Sample.females.Day.2
mali$total_asb_positive <- mali$females.ASB.positive 
mali$total_sampled <- mali$TOTAL.Sample.females.Day.2 
mali$total_catch <- mali$CDC.total.females.day.2
mali$days <- (mali$month-1)*30+mali$Day

par(las=1, mfrow=c(1,2), mar=c(8,4,4,1)+0.1)
plot(mali$days,
     mali$dyed_fraction,
     col=brewer.pal(7, "Set2")[factor(mali$Village)],
     cex=mali$total_sampled/mean(mali$total_sampled),
     pch=1,
     frame.plot=F,
     ylim=c(0,1),
     xlim=c(0,400),
     xlab = "Day",
     ylab = "Dyed fraction")
grid()
mali |>
  group_by(Village) |>
  summarise(dyed_fraction=weighted.mean(dyed_fraction,total_sampled),
            total_catch=sum(total_catch)) -> mali_grouped
par(las=2)
plot(mali_grouped$dyed_fraction,
     col=brewer.pal(7,"Set2")[factor(mali_grouped$Village)],
     cex=2*mali_grouped$total_catch/mean(mali_grouped$total_catch),
     frame.plot = F,
     ylim = c(0,1),
     ylab = "Dyed fraction",
     xaxt = "n",
     xlab = "")
axis(1, at=1:7, labels=mali_grouped$Village)
grid()

villages <- mali_grouped$Village
cluster_quantiles <- matrix(0,nrow = length(villages), ncol=2)
for (i in 1:length(villages)) {
  this_village <- villages[i]
  mali_filtered <- filter(mali, Village==this_village)
  replicates <- matrix(0,nrow=5000,ncol=nrow(mali_filtered))
  for (j in 1:5000) {
    replicates[j,] <- sample(x = mali_filtered$total_asb_positive/mali_filtered$total_sampled, 
                             replace = TRUE, 
                             size = nrow(mali_filtered), 
                             prob = mali_filtered$total_sampled)
  }
  cluster_means <- apply(replicates, MARGIN = 1, FUN = mean)
  cluster_quantiles[i,] <- quantile(cluster_means, probs = c(0.025,0.975))
}
colnames(cluster_quantiles) <- c("feed_rate_lower", "feed_rate_upper")
mali_grouped <- cbind(mali_grouped, cluster_quantiles)

arrows(1:7+0.15,
       mali_grouped[,4],
       1:7+0.15,
       mali_grouped[,5],
       code=0,
       lwd=2)

# cluster level abundance 
control_villages <- mali_grouped$Village
atsb_villages <- c("Tiko", "Kignele", "Krekrelo", "Farabale", "Trekrou", "Sirakele", "Sambadani")
villages <- c(control_villages, atsb_villages)
par(las=1, mfcol=c(7,2), mar = c(1.1, 1.1, 1.1, 1.1))
for (i in 1:14) {
  this_village <- villages[i]
  cdc_2016 |>
    filter(Vilage==this_village) -> this_cdc_2016
  cdc_2017 |>
    filter(Vilage==this_village) -> this_cdc_2017
  plot(this_cdc_2016$days,
       this_cdc_2016$total_catch,
       frame.plot = F,
       pch=20,
       xlab="Days",
       ylab="Catch",
       col=factor(this_cdc_2016$Location),
       xlim=c(0,800),
       ylim=c(0,max(c(this_cdc_2017$total_catch, this_cdc_2016$total_catch))))
  points(this_cdc_2017$days,
         this_cdc_2017$total_catch,
         pch=20)
  grid()
  title(this_village)
}

cdc_2017 |>
  group_by(days, Vilage, Experimental.or.control, Month) |>
  summarise(total_catch=sum(total_catch)) -> cdc_2017_grouped
cdc_2017_grouped$days <- cdc_2017_grouped$days-365
cdc_2016 |>
  group_by(year, days, Vilage, Location) |>
  summarise(total_catch=sum(total_catch)) |>
  filter(year == 16, Location=="inside", days > 90) -> cdc_2016_grouped


# baseline asb feeding in atsb clusters
library(DescTools)
atsb_villages <- c("Tiko", "Kignele", "Krekrelo", "Farabale", 
                   "Trekrou", "Sirakele", "Sambadani")

asb_feeding <- data_frame(
  village = c("Sirakele", "Farabale", "Tiko", "Sambadani", "Krekrelo", 
              "Trekrou", "Kignele"),
  date = c("06-10-2016", "08-10-2016", "19-08-2016", "26-08-2016", "09-09-2016", "20-09-2016", NA),
  total_catch = c(192, 87, 329, 216, 151, 181, NA),
  total_asb_positive = c(25+33, 33, 114, 69, 70, 69, NA)
)
asb_feeding$dyed_fraction <- asb_feeding$total_asb_positive/asb_feeding$total_catch

auc_2016 <- c()
auc_2017 <- c()
atsb_villages <- asb_feeding$village
control_villages <- mali_grouped$Village
villages <- c(atsb_villages, control_villages)
for (i in 1:length(villages)) {
  cdc_2016 |>
    group_by(year, Month, Vilage, Location) |>
    summarise(total_catch=sum(tot.f)) |>
    filter(Vilage==villages[i] & Location != "periphery" & year == 16) -> this_cdc_2016
  cdc_2017 |>
    group_by(Month, Vilage) |>
    summarise(total_catch=sum(tot.f)) |>
    filter(Vilage==villages[i]) -> this_cdc_2017
  auc_2016[i] <- AUC(x = this_cdc_2016$Month,
                     y = this_cdc_2016$total_catch,
                     from = 4,
                     to = 12)
  auc_2017[i] <- AUC(x = this_cdc_2017$Month,
                     y = this_cdc_2017$total_catch,
                     from = 4,
                     to = 12)
}

catch <- data_frame(village = villages,
                    auc_2016 = auc_2016,
                    auc_2017 = auc_2017,
                    ATSB = c(rep(1,7), rep(0,7)),
                    reduction = (auc_2016-auc_2017)/auc_2016)
catch$dyed_fraction <- c(asb_feeding$dyed_fraction, rep(NA,7))
catch$asb_sample <- c(asb_feeding$total_catch, rep(NA,7))

par(las=1)
plot(catch$dyed_fraction[1:6],
     catch$reduction[1:6],
     pch=1,
     col=factor(catch$village)[1:6],
     frame.plot = F,
     xlab = "2016 dyed fraction",
     ylab = "Reduction relative to 2016")
grid()
fit <- lm(catch$reduction[1:6] ~ catch$dyed_fraction[1:6])
summary(fit)
lines(catch$dyed_fraction[1:6],
      predict(fit))

cdc_2017_grouped |> 
  filter(Experimental.or.control == "Exp.") |>
  group_by(Month) |>
  summarise(total_catch=mean(total_catch)) -> atsb_2017
auc_control <- AUC(atsb_2017$Month,
                   atsb_2017$total_catch)
auc_2017_month <- c()
for (i in 1:length(villages)) {
  cdc_2017_grouped |>
    filter(Vilage==villages[i]) -> this_cdc_2017
  auc_2017_month[i] <- AUC(x = this_cdc_2017$Month,
                     y = this_cdc_2017$total_catch)
}
villages
catch$change_vs_control <- c(auc_2017_month[1:6]/auc_control, rep(NA, 8))

cdc_2017 |> 
  group_by(Month, Experimental.or.control) |>
  summarise(total_catch=sum(tot.f)/7) |>
  filter(Experimental.or.control == "Con.") -> control_2017
plot(control_2017$Month,
     control_2017$total_catch,
     frame.plot = F, 
     xlab = "Month",
     ylab = "Count",
     type = "l",
     lwd=2,
     ylim = c(0,600))
auc_control <- AUC(control_2017$Month,
                   control_2017$total_catch)
auc_village <- c()
for (i in 1:7) {
  this_village <- asb_feeding$village[i]
  cdc_2017 |>
    group_by(Month, Vilage) |>
    summarise(total_catch=sum(tot.f)) |>
    filter(Vilage == this_village) -> village_2017
  lines(village_2017$Month,
        village_2017$total_catch,
        lwd =2,
        col = brewer.pal(7,"Set1")[i])
  auc_village[i] <- AUC(village_2017$Month,
                        village_2017$total_catch)
}
grid()
reduction <- (auc_control - auc_village)/auc_control
plot(catch$dyed_fraction[1:6],
     reduction[1:6],
     col=brewer.pal(7, "Set1"),
     frame.plot = F,
     cex= catch$asb_sample[1:6]/mean(catch$asb_sample[1:6]),
     ylab = "Reduction relative to control",
     xlab = "Dyed fraction")
grid()
fit <- lm(reduction[1:6] ~ catch$dyed_fraction[1:6])
summary(fit)
lines(catch$dyed_fraction[1:6],
      predict(fit))


catch
