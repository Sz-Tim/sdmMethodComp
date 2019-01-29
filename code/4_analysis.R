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
walk(dir("code", "fn", full.names=T), source)
lam.df <- readRDS(here("vs", sp, "lam_df.rds"))

f_P <- list.files(here("out", sp), pattern="P_", full.names=T)
i_P <- extract_SDM_details(f_P, "P")
out_P <- map2(f_P, i_P, ~readRDS(.x) %>% mutate(SDM=.y[1], issue=.y[2])) %>%
  do.call(bind_rows, .) %>%
  full_join(dplyr::select(lam.df, 
                          -one_of("x", "y", "x_y", "inbd", "lat", "lon", "id")), 
            ., by="id.in") %>%
  mutate(fate_S=case_when(Surv.S>0 & prP>=0.5 ~ "S:1 P:1",
                          Surv.S==0 & prP>=0.5 ~ "S:0 P:1",
                          Surv.S>0 & prP<0.5 ~ "S:1 P:0",
                          Surv.S==0 & prP<0.5 ~ "S:0 P:0"),
         fate_lam=case_when(lambda>=1 & prP>=0.5 ~ "S:1 P:1",
                            lambda<1 & prP>=0.5 ~ "S:0 P:1",
                            lambda>=1 & prP<0.5 ~ "S:1 P:0",
                            lambda<1 & prP<0.5 ~ "S:0 P:0"))

out_LL <- out_P %>% 
  mutate(logP.N=dbinom(Surv.S>0, 1, prP, log=T),
         logP.N=ifelse(is.infinite(logP.N), -10, logP.N),
         logP.lam=dbinom(lambda>1, 1, prP, log=T),
         logP.lam=ifelse(is.infinite(logP.lam), -10, logP.lam)) %>%
  group_by(SDM, issue) %>%
  summarise(LL.N=sum(logP.N),
            LL.lam=sum(logP.lam)) %>%
  gather(Boundary, LogLik, 3:4)
out_LL$Boundary <- factor(out_LL$Boundary, labels=c("lambda > 1", "N > 0"))

f_TSS <- list.files(here("out", sp), pattern="TSS_", full.names=T)
i_TSS <- extract_SDM_details(f_TSS, "TSS")
out_TSS <- rbind(map(f_TSS, ~readRDS(.x)[[1]]) %>%
                   map2_dfr(., i_TSS, ~tibble(TSS=.x, Boundary="N > 0", 
                                              SDM=.y[1], issue=.y[2])), 
                 map(f_TSS, ~readRDS(.x)[[2]]) %>%
                   map2_dfr(., i_TSS, ~tibble(TSS=.x, Boundary="lambda > 1", 
                                              SDM=.y[1], issue=.y[2])))

if(overwrite) {
  write_csv(out_P, here("out", sp, "out_P.csv"))
  write_csv(out_TSS, here("out", sp, "out_TSS.csv"))
  write_csv(out_LL, here("out", sp, "out_LL.csv"))
}

