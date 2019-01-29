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
sp <- c("barberry", "garlic_mustard")[1]
overwrite <- TRUE
n_core_iss <- 4  # number of issues to run in parallel
n_core_obs <- 6  # number of simulations to run in parallel for each issue
n_sim <- 1 # number of simulations per sample (mechanistic only)

# load workspace
pkgs <- c("gbPopMod", "tidyverse", "magrittr", "MuMIn", "here", "doSNOW")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
sp_i <- read.csv("data/species_3km.csv") %>% filter(Name==sp)
vs.dir <- paste0("vs/", sp_i$Num, "/")
p <- readRDS(here(vs.dir, "p.rds"))
N_init <- readRDS(here(vs.dir, "N_init.rds"))
sdd.pr <- readRDS(here(vs.dir, "sdd.rds"))
sdd.ji <- readRDS(here(vs.dir, "sdd_ji.rds"))
p.ji <- readRDS(here(vs.dir, "p_ji.rds"))
lam.df <- readRDS(here(vs.dir, "lam_df.rds"))

clim_X_correct <- readRDS(here(vs.dir, "clim_X.rds"))
L <- map(dir(here(vs.dir), "env_", full.names=T), readRDS) %>%
  setNames(str_remove(dir(here(vs.dir), "env_"), ".rds"))
n.cell <- nrow(L$env_in); n.grid <- nrow(L$env_rct)
clim_X_wrong <- paste0("bio10_", c(1, 5, 12, "prMay")) %>%
  magrittr::extract(! . %in% clim_X_correct)
issue_i <- read.csv("data/issues.csv", stringsAsFactors=F)
vars_correct <- rep(0, length(clim_X_correct)*2+2)
names(vars_correct) <- c("(Intercept)", "size", 
                         paste0(rep(clim_X_correct, each=2), c("", "_sq")))
vars_wrong <- rep(0, length(clim_X_wrong)*2+2)
names(vars_wrong) <- c("(Intercept)", "size", 
                       paste0(rep(clim_X_wrong, each=2), c("", "_sq")))
out.dir <- paste0("out/", sp_i$Num)
if(!dir.exists(here(out.dir))) dir.create(here(out.dir), recursive=T)


p.c <- makeCluster(n_core_iss); registerDoSNOW(p.c)
foreach(i=seq_along(issue_i$Issue), .packages=c("dismo", pkgs)) %dopar% {
  # load issues
  issue <- issue_i$Issue[i]
  samp_iss <- issue_i$Sampling[i]
  mod_iss <- issue_i$Modeling[i]
  
  # set variables & model formulas
  v <- m <- list(CA=NULL, IPM=NULL)
  if(mod_iss=="wrongCov") {
    vars <- vars_wrong
    v.size <- grep("size", names(vars))
    v.i <- grep("bio10", names(vars))
  } else {
    vars <- vars_correct
    v.size <- grep("size", names(vars))
    v.i <- seq_along(vars)[-c(1,v.size)]
  }
  v$Mx <- grep("_sq|_cu", names(vars)[v.i], invert=T, value=T)
  v$CA <- vars[c(1, v.i)]
  m$CA <- paste(c(names(v$CA)[-1], "(1|yr)"), collapse=" + ")
  v$IPM <- vars[c(1, v.size, v.i)]
  m$IPM <- paste(names(v$IPM)[-1], collapse=" + ")
  n_z <- rep(list(1+length(v.size)), 4)  # adds intercept
  names(n_z) <- c("s", "g", "fl", "seed")
  if(!is.null(p$germ_x)) {
    n_x <- rep(list(length(v.i)), length(n_z)+1)
    names(n_x) <- c(names(n_z), "germ")
    n_x$germ <- n_x$germ + 1
  } else {
    n_x <- rep(list(length(v.i)), length(n_z))
    names(n_x) <- names(n_z)
  }
  
  # fit MaxEnt
  if("rJava" %in% rownames(installed.packages()) && 
     issue %in% issue_i$Issue[c(1:4,8)]) {
    P_MxE <- fit_MxE(sp_i$Num, issue, samp_iss, lam.df, v$Mx)
    if(overwrite) {
      saveRDS(P_MxE$diag, here(out.dir, paste0("Diag_MxE_", issue, ".rds")))
      saveRDS(P_MxE$P_MxE, here(out.dir, paste0("P_MxE_", issue, ".rds")))
      saveRDS(P_MxE$TSS_MxE, here(out.dir, paste0("TSS_MxE_", issue, ".rds")))
    }
  }
  # fit CA-demographic
  P_CA <- fit_CA(sp, sp_i, samp_iss, mod_iss, p,
                 L$env_rct, L$env_rct_unscaled, lam.df, v$CA, m$CA,
                 N_init, sdd.pr, sdd.ji, p.ji, n_sim, n_core_obs)
  if(overwrite) {
    saveRDS(P_CA$diag_CAd, here(out.dir, paste0("Diag_CAd_", issue, ".rds")))
    saveRDS(P_CA$P_CAd, here(out.dir, paste0("P_CAd_", issue, ".rds")))
    saveRDS(P_CA$TSS_CAd, here(out.dir, paste0("TSS_CAd_", issue, ".rds")))
  }
  # fit IPM, CA-individual
  P_IPM <- fit_IPM(sp, sp_i, samp_iss, mod_iss, p,
                   L$env_rct_unscaled, lam.df, v$IPM, m$IPM, n_z, n_x,
                   N_init, sdd.ji, p.ji, n_sim, n_core_obs)
  if(overwrite) {
    saveRDS(P_IPM$diag, here(out.dir, paste0("Diag_IPM_", issue, ".rds")))
    saveRDS(P_IPM$P_CAi, here(out.dir, paste0("P_CAi_", issue, ".rds")))
    saveRDS(P_IPM$P_IPM, here(out.dir, paste0("P_IPM_", issue, ".rds")))
    saveRDS(P_IPM$TSS_IPM, here(out.dir, paste0("TSS_IPM_", issue, ".rds")))
    saveRDS(P_IPM$TSS_CAi, here(out.dir, paste0("TSS_CAi_", issue, ".rds")))
  }
}
stopCluster(p.c)



