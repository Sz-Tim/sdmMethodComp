# 4: Analyze fit models
# Comparison of SDM approaches
# Tim Szewczyk

# This script aggregates and analyzes the output from `3_fitModels.R`.

########
## Setup
########
# file specifications
sp <- "sp1"
overwrite <- TRUE

# load workspace
pkgs <- c("tidyverse", "magrittr", "stringr")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
lam.df <- here(readRDS(paste0("out/", sp, "_lam_df.rds")))
f_IPM <- list.files(here("out"), pattern=paste0(sp, "_P_IPM"), full.names=T)
i_IPM <- extract_SDM_issues(f_IPM, "IPM")
out_IPM <- map2(f_IPM, i_IPM, ~readRDS(here(.x)) %>% 
                  mutate(SDM="IPM", s_Iss=.y[1], m_Iss=.y[2])) %>%
  do.call("rbind", .) %>%
  full_join(select(lam.df, -c(10:15,17)), ., by="id.inbd")
out_IPM$outcome <- NA
out_IPM$outcome[out_IPM$Surv.S>0 & out_IPM$Surv.S.f>0.5] <- "true1pred1"
out_IPM$outcome[out_IPM$Surv.S==0 & out_IPM$Surv.S.f<0.5] <- "true0pred0"
out_IPM$outcome[out_IPM$Surv.S>0 & out_IPM$Surv.S.f<0.5] <- "true1pred0"
out_IPM$outcome[out_IPM$Surv.S==0 & out_IPM$Surv.S.f>0.5] <- "true0pred1"


if(overwrite) {
  write.csv(out_IPM, here(paste0("out/", sp, "_IPM.csv")))
}

