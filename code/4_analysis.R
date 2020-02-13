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

# Aggregate predictions
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
out_P$B.f[out_P$issue=="noSB"] <- NA

# Compare issues to none
out_P <- full_join(out_P, 
                   filter(out_P, issue=="none") %>% select(-c(1:79, 92)), 
                   by=c("id", "SDM"), suffix=c("", ".none")) %>%
  mutate(sp=sp)

# Calculate log likelihood
out_LL <- out_P %>% 
  mutate(logP.N=dbinom(Surv.S>0, 1, prP, log=T),
         logP.N=ifelse(is.infinite(logP.N), -10, logP.N),
         logP.lam=dbinom(lambda>1, 1, prP, log=T),
         logP.lam=ifelse(is.infinite(logP.lam), -10, logP.lam)) %>%
  group_by(SDM, issue) %>%
  summarise(LL.N=sum(logP.N),
            LL.lam=sum(logP.lam)) %>%
  gather(Boundary, LogLik, 3:4) %>%
  mutate(sp=sp)
out_LL$Boundary <- factor(out_LL$Boundary, labels=c("lambda > 1", "N > 0"))

# Aggregate TSS
f_TSS <- list.files(here("out", sp), pattern="TSS_", full.names=T)
i_TSS <- extract_SDM_details(f_TSS, "TSS")
out_TSS <- rbind(map(f_TSS, ~readRDS(.)[[1]]) %>%
                   map2_dfr(., i_TSS, 
                            ~tibble(TSS=map_dbl(.x, ~.$TSS), 
                                    sensitivity=map_dbl(.x, ~.$sensitivity), 
                                    specificity=map_dbl(.x, ~.$specificity), 
                                    Boundary="N > 0", 
                                    SDM=.y[1], issue=.y[2])), 
                 map(f_TSS, ~readRDS(.)[[2]]) %>%
                   map2_dfr(., i_TSS, 
                            ~tibble(TSS=map_dbl(.x, ~.$TSS), 
                                    sensitivity=map_dbl(.x, ~.$sensitivity), 
                                    specificity=map_dbl(.x, ~.$specificity), 
                                    Boundary="lambda > 1", 
                                    SDM=.y[1], issue=.y[2]))) %>%
  mutate(sp=sp)

# Calculate MaxEnt diagnostics
f_MxE <- list.files(here("out", sp), pattern="Diag_MxE", full.names=T)
i_MxE <- extract_SDM_details(f_MxE, "Diag")
diagRaw_MxE <- map(f_MxE, readRDS)
AUC_MxE <- map2_dfr(diagRaw_MxE, i_MxE,
                     ~tibble(AUC=c(map_dbl(.x[[1]], ~.@auc), 
                                   map_dbl(.x[[2]], ~.@auc)),
                             Boundary=rep(c("N > 0", "lambda > 1"), each=100),
                             issue=.y[2])) %>%
  mutate(sp=sp)
ROC_MxE <- rbind(map2_dfr(diagRaw_MxE, i_MxE,
                          ~tibble(TPR=unlist(map(.x[[1]], ~.@TPR)),
                                  FPR=unlist(map(.x[[1]], ~.@FPR)),
                                  obs=rep(1:100, times=map_int(.x[[1]], 
                                                               ~length(.@TPR))),
                                  issue=.y[2],
                                  Boundary="N > 0")),
                 map2_dfr(diagRaw_MxE, i_MxE,
                          ~tibble(TPR=unlist(map(.x[[2]], ~.@TPR)),
                                  FPR=unlist(map(.x[[2]], ~.@FPR)),
                                  obs=rep(1:100, times=map_int(.x[[2]], 
                                                               ~length(.@TPR))),
                                  issue=.y[2],
                                  Boundary="lambda > 1"))) %>%
  mutate(sp=sp)


if(overwrite) {
  write_csv(out_P, here("out", sp, "out_P.csv"))
  write_csv(out_LL, here("out", sp, "out_LL.csv"))
  write_csv(out_TSS, here("out", sp, "out_TSS.csv"))
  write_csv(AUC_MxE, here("out", sp, "out_AUC.csv"))
  write_csv(ROC_MxE, here("out", sp, "out_ROC.csv"))
}

