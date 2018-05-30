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
n_core_iss <- 3  # number of issues to run in parallel
n_core_sim <- 8  # number of simulations to run in parallel for each issue
n_sim <- 16 # number of simulations per sample (mechanistic only)
vars <- c("(Intercept)"=0, 
          "temp"=0, "temp2"=0, "prec"=0, "prec2"=0, 
          "pOpn"=0, "pOth"=0, "pDec"=0, #"pEvg"=0, "pMxd"=0,
          "size"=0, "size2"=0)#, "size3"=0)

# load workspace
pkgs <- c("dismo", "gbPopMod", "tidyverse", "magrittr", "MuMIn", "here", "doSNOW")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), ~source(here(.)))
p <- readRDS(here(paste0("out/", sp, "_p.rds")))
N_init <- readRDS(here(paste0("out/", sp, "_N_init.rds")))
S <- readRDS(here(paste0("out/", sp, "_S.rds")))
U <- readRDS(here(paste0("out/", sp, "_U.rds")))
sdd.pr <- readRDS(here(paste0("out/", sp, "_sdd.rds")))
lam.df <- readRDS(here(paste0("out/", sp, "_lam_df.rds")))
env.in <- readRDS(here(paste0("out/", sp, "_env_in.rds")))
env.rct <- readRDS(here(paste0("out/", sp, "_env_rct.rds")))
env.rct.unsc <- readRDS(here(paste0("out/", sp, "_env_rct_unscaled.rds")))
n.cell <- nrow(env.in); n.grid <- nrow(env.rct)
issue_i <- read.csv(here("data/issues.csv"), stringsAsFactors=F)
S_p <- as.matrix(lam.df[lam.df$Surv.S>0, 15:14]) # true presences
S_a <- as.matrix(lam.df[lam.df$Surv.S==0, 15:14]) # true absences


p.c <- makeCluster(n_core_iss); registerDoSNOW(p.c)
foreach(i=seq_along(issue_i$Issue)) %dopar% {
  # load issues
  issue <- issue_i$Issue[i]
  samp_iss <- issue_i$Sampling[i]
  mod_iss <- issue_i$Modeling[i]
  
  # set model details
  v <- m <- n <- list(CA=NULL, IPM=NULL)
  v.size <- grep("size", names(vars))
  if(mod_iss=="clim") {
    v.i <- grep("temp|prec", names(vars))
  } else if(mod_iss=="lc") {
    v.i <- grep("^p(?!.*rec)", names(vars), perl=TRUE)
  } else {
    v.i <- seq_along(vars)[-c(1,v.size)]
  }
  v$CA <- vars[c(1, v.i)]
  m$CA <- paste(names(v$CA)[-1], collapse=" + ")
  v$IPM <- vars[c(1, v.size, v.i)]
  m$IPM <- paste(names(v$IPM)[-1], collapse=" + ")
  n$IPM$z <- rep(list(1+length(v.size)), 4)  # adds intercept
  n$IPM$x <- rep(list(length(v.i)), 4)
  names(n$IPM$z) <- names(n$IPM$x) <- c("s", "g", "fl", "seed")
  
  # fit MaxEnt
  P_Mx <- fit_Mx(sp, samp_iss, lam.df, vars, v.i, S_p, S_a)
  if(overwrite) {
    saveRDS(P_Mx$diag, here(paste0("out/", sp, "_Diag_Mx_", issue, ".rds")))
    saveRDS(P_Mx$pred, here(paste0("out/", sp, "_P_Mx_", issue, ".rds")))
  }
  # fit CA-demographic, CA-lambda
  P_CA <- fit_CA(sp, samp_iss, mod_iss, p, env.rct, env.rct.unsc, lam.df, vars,
                 v.i, v$CA, m$CA, N_init, sdd.pr, n.cell, n.grid, n_sim, n_core_sim)
  if(overwrite) {
    saveRDS(P_CA$diag_CAd, here(paste0("out/", sp, "_Diag_CAd_", issue, ".rds")))
    saveRDS(P_CA$diag_CAl, here(paste0("out/", sp, "_Diag_CAl_", issue, ".rds")))
    saveRDS(P_CA$P_CAd, here(paste0("out/", sp, "_P_CAd_", issue, ".rds")))
    saveRDS(P_CA$P_CAl, here(paste0("out/", sp, "_P_CAl_", issue, ".rds")))
  }
  # fit IPM
  P_IPM <- fit_IPM(sp, samp_iss, mod_iss, p, env.rct.unsc, lam.df, vars, v.i,
                   v$IPM, m$IPM, n$IPM$x, n$IPM$z, N_init, sdd.pr, n.cell,
                   n_sim, n_core_sim)
  if(overwrite) {
    saveRDS(P_IPM$diag, here(paste0("out/", sp, "_Diag_IPM_", issue, ".rds")))
    saveRDS(P_IPM$P_CAi, here(paste0("out/", sp, "_P_CAi_", issue, ".rds")))
    saveRDS(P_IPM$P_IPM, here(paste0("out/", sp, "_P_IPM_", issue, ".rds")))
  }
}
stopCluster(p.c)



