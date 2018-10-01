# 3: Fit models from samples
# Comparison of SDM approaches
# Tim Szewczyk

# This script uses the samples generated in `2_sampleDistribution.R` to fit a
# correlative MaxEnt SDM, a mechanistic CA-style SDM, and a mechanistic 
# IPM-style SDM. Model selection is performed for each model

########
## Setup
########
# file specifications
sp <- "sp1"
overwrite <- TRUE
n_core_iss <- 2  # number of issues to run in parallel
n_core_sim <- 2  # number of simulations to run in parallel for each issue
n_sim <- 2 # number of simulations per sample (mechanistic only)

# load workspace
pkgs <- c("dismo", "gbPopMod", "tidyverse", "magrittr", "MuMIn", "here", "doSNOW")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), ~source(here(.)))
p <- readRDS(here("vs", sp, "p.rds"))
N_init <- readRDS(here("vs", sp, "N_init.rds"))
sdd.pr <- readRDS(here("vs", sp, "sdd.rds"))
sdd.ji <- readRDS(here("vs", sp, "sdd_ji.rds"))
p.ji <- readRDS(here("vs", sp, "p_ji.rds"))
lam.df <- readRDS(here("vs", sp, "lam_df.rds"))
env.in <- readRDS(here("vs", sp, "env_in.rds"))
env.rct <- readRDS(here("vs", sp, "env_rct.rds"))
env.rct.unsc <- readRDS(here("vs", sp, "env_rct_unscaled.rds"))
n.cell <- nrow(env.in); n.grid <- nrow(env.rct)
issue_i <- read.csv(here("data/issues.csv"), stringsAsFactors=F)
vars <- rep(0, ncol(env.in)-9+4)
names(vars) <- c("(Intercept)", "size", "size2", "size2",
                 names(env.in)[1:(ncol(env.in)-9)])
vars <- vars[-c(4, 11:12, 16:17)]
if(!dir.exists(here("out", sp))) dir.create(here("out", sp), recursive=T)


p.c <- makeCluster(n_core_iss); registerDoSNOW(p.c)
foreach(i=seq_along(issue_i$Issue), .packages=pkgs) %dopar% {
  # load issues
  issue <- issue_i$Issue[i]
  samp_iss <- issue_i$Sampling[i]
  mod_iss <- issue_i$Modeling[i]
  
  # set model details
  v <- m <- list(CA=NULL, IPM=NULL)
  v.size <- grep("size", names(vars))
  if(mod_iss=="clim") {
    v.i <- grep("bio10", names(vars))
  } else if(mod_iss=="lc") {
    v.i <- grep("bio10", names(vars), invert=T)[-1]
  } else {
    v.i <- seq_along(vars)[-c(1,v.size)]
  }
  v$Mx <- grep("_sq|_cu", names(vars)[v.i], invert=T, value=T)
  v$CA <- vars[c(1, v.i)]
  m$CA <- paste(c(names(v$CA)[-1], "(1|yr)"), collapse=" + ")
  v$IPM <- vars[c(1, v.size, v.i)]
  m$IPM <- paste(names(v$IPM)[-1], collapse=" + ")
  n_z <- rep(list(1+length(v.size)), 4)  # adds intercept
  n_x <- rep(list(length(v.i)), 4)
  names(n_z) <- names(n_x) <- c("s", "g", "fl", "seed")
  
  # fit MaxEnt
  # if(issue %in% issue_i$Issue[c(1:4,8:9)]) {
  #   P_MxE <- fit_MxE(sp, issue, samp_iss, lam.df, v$Mx)
  #   if(overwrite) {
  #     saveRDS(P_MxE$diag, here("out", sp, paste0("Diag_MxE_", issue, ".rds")))
  #     saveRDS(P_MxE$P_MxE, here("out", sp, paste0("P_MxE_", issue, ".rds")))
  #   }
  # }
  # fit CA-demographic
  P_CA <- fit_CA(sp, samp_iss, mod_iss, p, env.rct, env.rct.unsc, lam.df, 
                 v$CA, m$CA, N_init, sdd.pr, sdd.ji, p.ji, n_sim, n_core_obs)
  if(overwrite) {
    saveRDS(P_CA$diag_CAd, here("out", sp, paste0("Diag_CAd_", issue, ".rds")))
    saveRDS(P_CA$P_CAd, here("out", sp, paste0("P_CAd_", issue, ".rds")))
    saveRDS(P_CA$P_CAl, here("out", sp, paste0("P_CAl_", issue, ".rds")))
  }
  # fit IPM, CA-individual
  P_IPM <- fit_IPM(sp, samp_iss, mod_iss, p, env.rct.unsc, lam.df, v$IPM, m$IPM, 
                   n_z, n_x, N_init, sdd.ji, p.ji, n_sim, n_core_obs)
  if(overwrite) {
    saveRDS(P_IPM$diag, here("out", sp, paste0("Diag_IPM_", issue, ".rds")))
    saveRDS(P_IPM$P_CAi, here("out", sp, paste0("P_CAi_", issue, ".rds")))
    saveRDS(P_IPM$P_IPM, here("out", sp, paste0("P_IPM_", issue, ".rds")))
  }
}
stopCluster(p.c)



