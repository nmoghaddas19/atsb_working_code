# Here I'm going to use some of the catch data from Mali and Zambia to 
# estimate phi_indoors which is a parameter describing the fraction of 
# mosquito bites taken indoors. Malariasimulation uses phi_indoors as a
# parameter to modulate the efficacy of certain interventions such as 
# bednets and IRS. Phi_indoors is ideally calculated using human landing
# catch data but can also be done using trap data such as UV light traps
# if we assume that every mosquito caught in a trap was making an attempt to 
# blood feed. The method for calculating phi_indoors is described in 
# Sherrard-Smith et al., (2019). It uses both the hourly catch data and the 
# distribution of the proportion of people outdoors in a given hour. The later
# I believe was not recorded in the ASB/ATSB trials but there may be published 
# results for Mali and Zambia on the literature. 

# load zambia data
zambia <- read.csv("atsb_working_code/ZAM.ASB.Target Summaries.2021.10.01.clean.csv")

# there are columns in the data frame for whether each mosquito was caught
# indoors or outdoors
table(zambia$Indoor)
table(zambia$Outdoor)
24675/(24675+19845) # looks like ~ 55% of bites were taken indoors overall

# the data frame does not have information on the hour that mosquitoes were 
# caught so this is as far as we can go 
table(zambia$An..gambiae)
table(zambia$An..funestus)

# load mali data
hlc_apr_jul_2017 <- read.csv("atsb_working_code/HLC data processed/April_May_June_July-Table 1.csv")[,1:16]
hlc_aug_sep_2017 <- read.csv("atsb_working_code/HLC data processed/August_September-Table 1.csv")[,1:16]
hlc_oct_dec_2017 <- read.csv("atsb_working_code/HLC data processed/Oct_Nvber_Dcber-Table 1.csv")[,1:16]
hlc <- rbind(hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)
rm(hlc_apr_jul_2017, hlc_aug_sep_2017, hlc_oct_dec_2017)

# there is a column in the dataframe for position
table(hlc$Position)
# seems like they couldn't decide between outside and Outdoor

