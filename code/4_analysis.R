# 4: Analyze fit models
# Comparison of SDM approaches
# Tim Szewczyk

# This script aggregates and analyzes the output from `3_fitModels.R`.

########
## Setup
########
# file specifications
sp <- "sp1"

# load workspace
pkgs <- c("tidyverse", "magrittr", "stringr")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
source("code/fn_IPM.R"); source("code/fn_aux.R"); source("code/fn_sim.R")
lam.df <- readRDS(paste0("out/", sp, "_lam_df.rds"))
f_IPM <- list.files("out", pattern=paste0(sp, "_P_IPM"), full.names=T)
i_IPM <- extract_SDM_issues(f_IPM, "IPM")
out_IPM <- map2(f_IPM, i_IPM, ~readRDS(.x) %>% 
                  mutate(SDM="IPM", s_Iss=.y[1], m_Iss=.y[2])) %>%
  do.call("rbind", .) %>%
  full_join(select(lam.df, -c(10:15,17)), ., by="id.inbd")

write.csv(out_IPM, paste0("out/", sp, "_IPM.csv"))

