# This script will draw heavily on my ASB workflow to generate estimates 
# of the uncertainty in prevalence caused by not measuring feeding rates
# in trial clusters. I will do this by assuming the variability in dyed 
# fraction in control clusters is comparable to what would be seen in the
# trial clusters. 

# read in the dyed fraction data
zambia <- read.csv("~/Documents/GitHub/atsb_working_code/zambia_asb_data")
mali <- read.csv("~/Documents/GitHub/atsb_working_code/mali_asb_data")

# extract static cluster and species specific dyed fractions
