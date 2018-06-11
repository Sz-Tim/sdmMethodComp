# 4: Analyze fit models
# Comparison of SDM approaches
# Tim Szewczyk

# This script aggregates and analyzes the output from `3_fitModels.R`.

########
## Setup
########
# file specifications
sp <- "sp3"
overwrite <- TRUE

# load workspace
pkgs <- c("tidyverse", "magrittr", "stringr", "here")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
lam.df <- readRDS(here(paste0("out/", sp, "_lam_df.rds")))

f_P <- list.files(here("out"), pattern=paste0(sp, "_P_"), full.names=T)
i_P <- extract_SDM_details(f_P)
out <- map2(f_P, i_P, ~readRDS(.x) %>% mutate(SDM=.y[1], issue=.y[2])) %>%
  do.call(bind_rows, .) %>%
  full_join(dplyr::select(lam.df, -c(10:15,17)), ., by="id.inbd") %>%
  mutate(boundary=case_when(SDM %in% c("CAd", "CAi", "MxE", "MxL") ~ "Surv",
                            SDM == "IPM" ~ "lam"),
         outcome=case_when(boundary=="Surv" & Surv.S>0 & prP>=0.5 ~ "S:1 P:1",
                           boundary=="Surv" & Surv.S==0 & prP>=0.5 ~ "S:0 P:1",
                           boundary=="Surv" & Surv.S>0 & prP<0.5 ~ "S:1 P:0",
                           boundary=="Surv" & Surv.S==0 & prP<0.5 ~ "S:0 P:0",
                           boundary=="lam" & lambda>=1 & prP>=0.5 ~ "S:1 P:1",
                           boundary=="lam" & lambda<1 & prP>=0.5 ~ "S:0 P:1",
                           boundary=="lam" & lambda>=1 & prP<0.5 ~ "S:1 P:0",
                           boundary=="lam" & lambda<1 & prP<0.5 ~ "S:0 P:0"))

if(overwrite) {
  write.csv(out, here(paste0("out/", sp, "_out.csv")), row.names=F)
}

