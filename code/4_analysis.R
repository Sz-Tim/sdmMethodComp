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
lam.df <- readRDS(here("vs", sp, "lam_df.rds"))

f_P <- list.files(here("out", sp), pattern="P_", full.names=T)
i_P <- extract_SDM_details(f_P)
out <- map2(f_P, i_P, ~readRDS(.x) %>% mutate(SDM=.y[1], issue=.y[2])) %>%
  do.call(bind_rows, .) %>%
  full_join(dplyr::select(lam.df, -one_of("x", "y", "x_y", "inbd", "lat", "lon", "id")), 
            ., by="id.inbd") %>%
  mutate(fate_S=case_when(Surv.S>0 & prP>=0.5 ~ "S:1 P:1",
                          Surv.S==0 & prP>=0.5 ~ "S:0 P:1",
                          Surv.S>0 & prP<0.5 ~ "S:1 P:0",
                          Surv.S==0 & prP<0.5 ~ "S:0 P:0"),
         fate_lam=case_when(lambda>=1 & prP>=0.5 ~ "S:1 P:1",
                            lambda<1 & prP>=0.5 ~ "S:0 P:1",
                            lambda>=1 & prP<0.5 ~ "S:1 P:0",
                            lambda<1 & prP<0.5 ~ "S:0 P:0"))

if(overwrite) {
  write_csv(out, here("out", sp, "out.csv"))
}

