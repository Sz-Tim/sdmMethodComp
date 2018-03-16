# 5: Plot output
# Comparison of SDM approaches
# Tim Szewczyk

# This script plots the output from `3_fitModels.R` and `4_analysis.R`.

########
## Setup
########
# file specifications
sp <- "sp1"


# barplots showing % correct for presence/absence; fill=SDM; facet=issues
ggplot(out_IPM, aes(x=sign(Surv.S.f) == sign(Surv.S), fill=s_Iss)) + 
  geom_bar(position="dodge")