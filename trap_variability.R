

library(lme4)

# load in trap data here

cdc_2017_bytrap <- matrix(0, nrow=1260, ncol=6)

for (i in 1:nrow(cdc_2017)) {
  for (j in 1:10) {
    cdc_2017_bytrap[(i-1)*10+j,6] <- cdc_2017[i,j+7]
    cdc_2017_bytrap[(i-1)*10+j,5] <- j
    cdc_2017_bytrap[(i-1)*10+j,4] <- cdc_2017[i,4]
    cdc_2017_bytrap[(i-1)*10+j,3] <- cdc_2017[i,3]
    cdc_2017_bytrap[(i-1)*10+j,2] <- cdc_2017[i,2]
    cdc_2017_bytrap[(i-1)*10+j,1] <- cdc_2017[i,1]
  }
}
colnames(cdc_2017_bytrap) <- c("Month", "Date", "Vilage", "Experimental.or.control",
                               "Trap_number", "Count")
cdc_2017_bytrap <- data.frame(cdc_2017_bytrap)
cdc_2017_bytrap$Observation <- as.factor(1:1260)
cdc_2017_bytrap$Month <- as.numeric(cdc_2017_bytrap$Month)
cdc_2017_bytrap |>
  filter(Month > 6 & Month < 12) -> cdc_2017_bytrap

cdc_2017_bytrap$Month <- as.factor(cdc_2017_bytrap$Month)
cdc_2017_bytrap$Vilage <- as.factor(cdc_2017_bytrap$Vilage)
cdc_2017_bytrap$Trap_number <- as.factor(cdc_2017_bytrap$Trap_number)
cdc_2017_bytrap$Count <- as.numeric(cdc_2017_bytrap$Count)

glmer.nb(Count ~ as.factor(Experimental.or.control) +
           (1|Month) + (1|Vilage) + (1|Trap_number), data=cdc_2017_bytrap) |> 
  summary()

cdc_2017_bytrap |> 
  filter(Experimental.or.control=="Exp.") |>
  summarise(Mean=mean(Count))
 