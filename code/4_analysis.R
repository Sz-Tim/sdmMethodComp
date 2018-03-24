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
lam.df <- readRDS(here(paste0("out/", sp, "_lam_df.rds")))

f_CA <- list.files(here("out"), pattern=paste0(sp, "_P_CA"), full.names=T)
i_CA <- extract_SDM_issues(f_CA, "CA")
out_CA <- map2(f_CA, i_CA, ~readRDS(.x) %>% 
                  mutate(SDM="CA", s_Iss=.y[1], m_Iss=.y[2])) %>%
  do.call("rbind", .) %>%
  full_join(select(lam.df, -c(10:15,17)), ., by="id.inbd")
out_CA$outcome <- NA
out_CA$outcome[out_CA$Surv.S>0 & out_CA$Surv.S.f>0.5] <- "true1pred1"
out_CA$outcome[out_CA$Surv.S==0 & out_CA$Surv.S.f<0.5] <- "true0pred0"
out_CA$outcome[out_CA$Surv.S>0 & out_CA$Surv.S.f<0.5] <- "true1pred0"
out_CA$outcome[out_CA$Surv.S==0 & out_CA$Surv.S.f>0.5] <- "true0pred1"


f_IPM <- list.files(here("out"), pattern=paste0(sp, "_P_IPM"), full.names=T)
i_IPM <- extract_SDM_issues(f_IPM, "IPM")
out_IPM <- map2(f_IPM, i_IPM, ~readRDS(.x) %>% 
                  mutate(SDM="IPM", s_Iss=.y[1], m_Iss=.y[2])) %>%
  do.call("rbind", .) %>%
  full_join(select(lam.df, -c(10:15,17)), ., by="id.inbd")
out_IPM$outcome <- NA
out_IPM$outcome[out_IPM$Surv.S>0 & out_IPM$Surv.S.f>0.5] <- "true1pred1"
out_IPM$outcome[out_IPM$Surv.S==0 & out_IPM$Surv.S.f<0.5] <- "true0pred0"
out_IPM$outcome[out_IPM$Surv.S>0 & out_IPM$Surv.S.f<0.5] <- "true1pred0"
out_IPM$outcome[out_IPM$Surv.S==0 & out_IPM$Surv.S.f>0.5] <- "true0pred1"

IPM_only <- names(out_IPM)[!names(out_IPM) %in% names(out_CA)]
out <- bind_rows(out_IPM, out_CA)


if(overwrite) {
  write.csv(out, here(paste0("out/", sp, "out.csv")))
  write.csv(out_CA, here(paste0("out/", sp, "out_CA.csv")))
  write.csv(out_IPM, here(paste0("out/", sp, "out_IPM.csv")))
}

