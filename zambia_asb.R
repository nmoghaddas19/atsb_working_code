# The aim of this script is to examine the ASB trial data from Zambia to 
# parameterise malariasimulation and compare predictions to observed values 
# for mosquito counts (are there in this dataset?)

library(lubridate)
library(RColorBrewer)
library(lme4)
# read in data
zambia <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")

# looking at the data it seems like there is no ATSB arm. all clusters had either 
# 2 or 3 ASBs deployed per eligible structure throughout so I cant compare predicted
# counts to anything if i generate them

# i guess we can still have a look at the feeding rate. first we need to collapse
# the data frame

zambia |> 
  mutate(Month=month(parse_date_time(Day.of.Date.Trap, "dmy"))) |>
  group_by(Month, Cluster) |>
  summarise(ASB_fraction=sum(form.DyeFed)/n(), 
            Count = n()) -> zambia

d_asb <- matrix(0, nrow=3, ncol=3)
for (i in 3:5) {
  zambia |>
    filter(Month == i) -> t
  x <- matrix(rep(t$ASB_fraction, 5000), byrow = T, ncol = length(t$ASB_fraction))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,length(t$ASB_fraction),T))-mean(x)})
  d_asb[i-2,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$ASB_fraction)
}  
par(mfrow=c(1,1), las=1)
plot(zambia$Month,
     zambia$ASB_fraction,
     pch=20,
     frame.plot = F, ylim=c(0,0.35))
polygon(c(3:5,5:3), c(d_asb[,1],rev(d_asb[,3])),
        border = FALSE, col = adjustcolor("black", alpha.f = 0.2))

d_asb <- matrix(0, nrow=3, ncol=3)
for (i in 3:5) {
  zambia |>
    filter(Month == i) -> t
  x <- matrix(rep(t$Count, 5000), byrow = T, ncol = length(t$Count))
  a <- apply(x, MARGIN = 1, FUN = function(x){mean(sample(x,length(t$Count),T))-mean(x)})
  d_asb[i-2,] <- quantile(a, c(0.025,0.5,0.975)) + mean(t$Count)
}  
plot(zambia$Month,
     zambia$Count,
     pch=20,
     frame.plot = F)
polygon(c(3:5,5:3), c(d_asb[,1],rev(d_asb[,3])),
        border = FALSE, col = adjustcolor("black", alpha.f = 0.2))

# by day
zambia <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")

zambia |> 
  filter(An..funestus==1) |>
  filter(Cluster == "asb84") |>
  mutate(Date=parse_date_time(Day.of.Date.Trap, "dmy")) |>
  group_by(Date) |>
  summarise(ASB_fraction=sum(form.DyeFed)/n(), 
            Count = n()) -> zambia
  
plot(zambia$Date,
     zambia$Count,
     pch=20,
     cex=zambia$ASB_fraction*3,
     frame.plot = F, xlab = "Date", ylab = "Count")
plot(zambia$Date,
     zambia$ASB_fraction,
     pch=20,
     frame.plot = F, xlab = "Date", ylab = "Dyed fraction",
     cex=zambia$Count/50,
     ylim=c(0,1))

# by day by cluster
zambia <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")
zambia |> 
  mutate(Date=parse_date_time(Day.of.Date.Trap, "dmy")) |>
  group_by(Date, Cluster) |>
  summarise(ASB_fraction=sum(form.DyeFed)/n(), 
            Count = n()) -> zambia
table(zambia$Cluster, zambia$Date)
plot(zambia$Date,
     zambia$Count,
     pch=20,
     cex=1.2,
     col=brewer.pal(10, "Set3")[as.factor(zambia$Cluster)],
     frame.plot = F, xlab = "Date", ylab = "Count")

plot(zambia$Date,
     zambia$ASB_fraction,
     pch=20,
     cex=1.2,
     col=brewer.pal(10, "Set3")[as.factor(zambia$Cluster)],
     frame.plot = F, xlab = "Date", ylab = "Dyed fraction")

# Functions ####
my.barplot <-
  function(tab, percent=F, ...) {
    x <- barplot(tab, beside = TRUE, legend = TRUE, col = grey((nrow(tab):1)/(0.5 + nrow(tab))),
                 ylim = c(0, max(tab) * 1.15), las = 2, ...)
    abline(h = 0)
    lab <- if(percent) sprintf("%.1f", tab) else as.character(tab)
    if(percent) lab <- paste0(lab, "%")
    text(x,tab,labels = lab, srt = 90, cex = 0.7, offset = 10,
         adj = -0.2)
    return(NULL)
  }

jensen.logit.adjust <- 
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

binom.glmm.test.binary.factor <-
  function(term, rest.of.formula, group.by, data, check.resid.plot = FALSE,
           add.odds.ratio = TRUE,
           feed.rate.pars = NULL, ...) { # extract levels
    levs <- levels(data[, term])
    # make formula and fit binomial GLMM
    form <- paste(rest.of.formula, term, sep = " + ")
    fit <- glmer(form, family = "binomial", data = data, control = glmerControl(optimizer = "bobyqa"))
    print(summary(fit))
    if(check.resid.plot) GLMMmisc::sim.residplot(fit)
    # get p-value for effect of term
    #p.val <- drop1(fit, test = "Chisq")[term, "Pr(Chi)"]
    #p.val.fmt <- if(p.val < 0.0005) "P < 0.001" else paste("P =", sprintf("%.3f", p.val)) 
    # get predictions for levels of term (only works if there are no fixed effects except term 
    fit0 <- update(fit, ~ . -1)
    est.tab.bias <-
      cbind(fixef(fit0),
            confint(fit0, method = "Wald")[paste0(term, levs), ])
    # adjust predictions for Jensen's inequality
    est.tab <-
      jensen.logit.adjust(p = plogis(est.tab.bias), V = sum(unlist(VarCorr(fit0))))
    # get effect estimate for non-reference level
    eff.tab <-
      c(fixef(fit)[paste0(term, levs[2])],
        confint(fit, method = "Wald")[paste0(term, levs[2]), ])
    or.tab <- sprintf("%.2f", exp(eff.tab))
    or.tab.fmt <- paste0(or.tab[1], " (", or.tab[2], ", ", or.tab[3], ")")
    # make plot with these results
    mf <- model.frame(fit)
    tab <- table(mf[, term], mf[, group.by], mf[, 1])
    proptab <- prop.table(tab, 1:2)[, , 2]
    propdf <- na.omit(data.frame(t(proptab)))
    propdf$n <- table(mf[, group.by])[as.character(propdf$Var1)]
    # make plot
    ymax <- max(c(propdf$Freq, c(est.tab)))
    stripchart(Freq ~ Var2, axes = TRUE, vertical = TRUE, method = "jitter", jitter = 0.07, 
               pch = 21, data = propdf, xlim = c(0.5, nrow(tab) + 0.5),
               ylim = c(0, 1.5 * ymax),
               cex = 2 * sqrt(propdf$n)/max(sqrt(propdf$n)), ...)
    if(add.odds.ratio) {
      lines(x = c(1, nrow(tab)), y = rep(1.2 * ymax, nrow(tab)), lwd = 0.7) 
      lines(x = c(1, 1), y = ymax * c(1.2, 1.18), lwd = 0.7)
      lines(x = c(nrow(tab), nrow(tab)), y = ymax * c(1.2, 1.18), lwd = 0.7) 
      text(x = mean(c(1, nrow(tab))), y = rep(1.2 * ymax, nrow(tab)),
           labels = p.val.fmt, pos = 3)
      text(x = mean(c(1, nrow(tab))), y = rep(1.3 * ymax, nrow(tab)),
           labels = paste("OR (95% CI) =", or.tab.fmt),
           pos = 3) }
    arrows(x0 = c(1, nrow(tab)) + 0.15, x1 = c(1, nrow(tab)) + 0.15,
           y0 = est.tab[, "2.5 %"], y1 = est.tab[, "97.5 %"],
           angle = 90, code = 3, length = 0.1)
    points(c(1, nrow(tab)) + 0.15, est.tab[, 1], pch = 23, bg = "grey")
    text(x = c(1, nrow(tab)) + 0.15, y = est.tab[, 1],
         labels = sprintf("%.3f", est.tab[, 1]),
         pos = 4)
    #return(fit)
  }
    
    
    
# stop ####

dat.all <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv",
                    stringsAsFactors = FALSE)
names(dat.all) <- tolower(names(dat.all))
dim(dat.all)
table(table(dat.all$form.sampleid.check))
dat.all <- dat.all[!duplicated(dat.all$form.sampleid.check), ]
table(table(dat.all$form.sampleid.check))
dat.all$mosquito_species[dat.all$an..gambiae %in% 1] <- "gambiae"
dat.all$mosquito_species[dat.all$an..funestus %in% 1] <- "funestus"
dat.all$mosquito_species <- factor(dat.all$mosquito_species)
dat.all$location <- factor(dat.all$indoor, 0:1, c("Outdoors", "Indoors"))
dat.all$anoph_sex <- factor(dat.all$male, 0:1, c("female", "male"))
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
# dat.all$form.bloodfed <- factor(dat.all$form.bloodfed)

sp <- "funestus"
dat <- droplevels(dat.all[dat.all$mosquito_species  %in% sp, ])
n.days <- as.vector(diff(range(dat$date.id))) + 1
hist(dat$date.id, breaks = n.days, freq = TRUE,
     main = paste("Daily numbers of female An.", sp, "collected over",
                  round(n.days/7, 1), "weeks\nTotal number collected:", nrow(dat)),
     col = "grey", xlab = "", las = 2)

my.barplot(tab = table(dat$location, dat$cluster),
           main = paste("An.", sp, "females collected by cluster\n",
                        nrow(dat),"females collected over",  nlevels(dat$cluster.day),
                        "collections"),
           ylab = "Total number collected")
my.barplot(tab = 100 * prop.table(table(dat$location, dat$cluster, dat$positive), 1:2)[, , 2],
           main = paste("An.", sp, "female dye positivity by cluster"),
           ylab = "Proportion (%) dye-positive",
           percent = TRUE)

#' Make plot showing change in positivity rate over time
dat.by.date.tab <- do.call("rbind", tapply(factor(dat$positive), dat$date.id, table))
dat.by.date <- as.data.frame(dat.by.date.tab)
dat.by.date$date <- as.Date(rownames(dat.by.date))
dat.by.date$n <- dat.by.date$`0` + dat.by.date$`1`
dat.by.date$prop.positive <- 100 * dat.by.date$`1` / dat.by.date$n
plot(prop.positive ~ date, data = dat.by.date, pch = 21, bg = "grey", xlab = "",
     cex = 2 * sqrt(n)/max(sqrt(n)), type = "b", ylab = "Proportion (%) dye-positive")
title(paste("Daily dye-positive proportion among female An.", sp))
cex.scale <- c(10, 20, 50, 100, 200, 500, 1000, 2000)
cex.scale <- cex.scale[cex.scale < max(dat.by.date$n)]
legend("topright", legend = paste0("n=", cex.scale), pch = 21, pt.bg = "grey",
       pt.cex = 2 * sqrt(cex.scale)/max(sqrt(dat.by.date$n)))



# Regression modelling
# dat |> 
#   filter(!form.bloodfed==3) -> dat
rest.of.formula <- "positive ~ (1|hh) + (1 | collection_date) + (1|cluster)"
        # "positive ~ (1 | cluster) + (1 | hh) + (1 | collection_date) + (1 | cluster.day)")
binom.glmm.test.binary.factor(term = "form.bloodfed", rest.of.formula = rest.of.formula,
                              group.by = "cluster", data = dat, ylab = "Proportion dye-positive",
                              main = paste("Proportion of dye-positive An.", sp, "\nby location and cluster"))

term <- "form.bloodfed"
rest.of.formula <- "positive ~ (1|hh) + (1 | collection_date) + (1|cluster)"
form <- paste(rest.of.formula, term, sep = " + ")                             
fit <- glmer(form, family = "binomial", data = dat, 
             control = glmerControl(optimizer = "bobyqa"))
print(summary(fit))

# get p-value for effect of term
p.val <- drop1(fit, scope=c("cluster","month_caught", "cluster:month_caught") , test = "Chisq")[term, "Pr(Chi)"]

# get predictions for levels of term (only works if there are no fixed effects except term 
fit0 <- update(fit, ~ . -1)
levs <- c("0", "1", "3")
est.tab.bias <-
  cbind(fixef(fit0),
        confint(fit0, method = "Wald")[paste0(term, levs), ])
# adjust predictions for Jensen's inequality
est.tab <-
  jensen.logit.adjust(p = plogis(est.tab.bias), V = sum(unlist(VarCorr(fit0))))                       
      

dat |> 
  filter(!is.na(dat$form.bloodfed)) |>
  group_by(cluster, form.bloodfed) |>
  summarise(ASB_fraction=sum(positive)/n(), 
            Count = n()) -> zambia
plot(zambia$form.bloodfed,
     zambia$ASB_fraction,
     pch=20, frame.plot=F)
arrows(x0 = c(0,1,3),
       y0 = est.tab[,2],
       x1 = c(0,1,3),
       y1 = est.tab[,3],
       col=4,
       code=3,
       angle=90,
       length=0.1)

library(lme4)
par(mfrow=c(1,1))
sugar_feeding <- read.csv("atsb_working_code/DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")
sugar_feeding <- sugar_feeding[-(64:65),] # final two rows are NAs so removing them   
par(las=1)
sugar_feeding$dyed_fraction <- sugar_feeding$females.ASB.positive/sugar_feeding$TOTAL.Sample.females.Day.2
sugar_feeding$days <- (sugar_feeding$month-1)*30+sugar_feeding$Day
plot.default(sugar_feeding$days, 
             sugar_feeding$dyed_fraction,
             cex=1.1, pch=20, frame.plot = F, xlab = "Day", ylab = "% bait fed",
             ylim = c(0,1), cex.axis=1.2, xlim = c(1,365))
grid()

# We are interested to know whether feeding rates change over time which would
# provide reason for including a time varying feeding rate in the final 
# simulations. 
# FIRST APPROACH: Fit a glmm with random effects on cluster and date. If the 
# random effect for date is greater than cluster, this provides evidence that 
# there is variation in feeding rates over time over and above the variation
# between clusters

# SECOND APPROACH: Fit a glmm with a fixed effect on date. If there is a 
# significant effect of date, this provides evidence that the feeding rate is
# changing over time


fit <-
  glmer(
    cbind(total_sampled, total_sampled - total_asb_positive) ~ 
      Village + (1|days),
    family = "binomial", data = sugar_feeding,
    control = glmerControl(optimizer = "bobyqa"))
summary(fit)

sp <- c("funestus")
dat <- droplevels(dat.all[dat.all$mosquito_species  %in% sp, ])
term <- "cluster"
rest.of.formula <- "positive ~ (1 | collection_date) + (1|hh)"
form <- paste(rest.of.formula, term, sep = " + ")  


fit <- glmer(positive ~ (1|collection_date) + (1|cluster),
             family = "binomial", data = dat.all, 
             control = glmerControl(optimizer = "bobyqa"))
print(summary(fit))

fit <- glmer(
  positive ~ (1|cluster) + (1|days),
  data = zambia, family = "binomial", control = glmerControl(optimizer = "bobyqa")
)
summary(fit)
