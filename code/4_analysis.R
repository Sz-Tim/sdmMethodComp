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
pkgs <- c("tidyverse", "magrittr", "stringr", "here")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
lam.df <- readRDS(here(paste0("out/", sp, "_lam_df.rds")))

f_CA <- list.files(here("out"), pattern=paste0(sp, "_P_CA"), full.names=T)
i_CA <- extract_SDM_issues(f_CA, "CA")
out_CA <- map2(f_CA, i_CA, ~readRDS(.x) %>% 
                  mutate(SDM="CA", s_Iss=.y[1], m_Iss=.y[2])) %>%
  do.call("rbind", .) %>%
  full_join(select(lam.df, -c(10:15,17)), ., by="id.inbd")
out_CA$outcome <- with(out_CA, case_when(
  Surv.S>0 & Surv.S.f>0.5 ~ "S:1 P:1",
  Surv.S==0 & Surv.S.f>0.5 ~ "S:0 P:1",
  Surv.S>0 & Surv.S.f<0.5 ~ "S:1 P:0",
  Surv.S==0 & Surv.S.f<0.5 ~ "S:0 P:0"
))


f_IPM <- list.files(here("out"), pattern=paste0(sp, "_P_IPM"), full.names=T)
i_IPM <- extract_SDM_issues(f_IPM, "IPM")
out_IPM <- map2(f_IPM, i_IPM, ~readRDS(.x) %>% 
                  mutate(SDM="IPM", s_Iss=.y[1], m_Iss=.y[2])) %>%
  do.call("rbind", .) %>%
  full_join(select(lam.df, -c(10:15,17)), ., by="id.inbd")
out_IPM$outcome <- with(out_IPM, case_when(
  Surv.S>0 & Surv.S.f>0.5 ~ "S:1 P:1",
  Surv.S==0 & Surv.S.f>0.5 ~ "S:0 P:1",
  Surv.S>0 & Surv.S.f<0.5 ~ "S:1 P:0",
  Surv.S==0 & Surv.S.f<0.5 ~ "S:0 P:0"
))

IPM_only <- names(out_IPM)[!names(out_IPM) %in% names(out_CA)]
out <- bind_rows(out_IPM, out_CA)


if(overwrite) {
  write.csv(out, here(paste0("out/", sp, "out.csv")))
  write.csv(out_CA, here(paste0("out/", sp, "out_CA.csv")))
  write.csv(out_IPM, here(paste0("out/", sp, "out_IPM.csv")))
}

