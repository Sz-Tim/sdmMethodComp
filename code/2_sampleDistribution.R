# 2: Sample from simulated distribution
# Comparison of SDM approaches
# Tim Szewczyk

# This script samples from the distribution simulated in `1_simulateSpecies.R` 
# and formats the samples for the purposes of fitting a correlative Maxent SDM,
# a mechanistic CA-style SDM, and a mechanistic IPM-style SDM. The issues we
# explore that are sampling-based are imposed here.

########
## Setup
########
# file specifications
sp <- c("barberry", "garlic_mustard")[1]
overwrite <- TRUE
samp.issues <- c("none", "noise", "sampBias", "nonEq")

# load workspace
pkgs <- c("tidyverse", "magrittr", "here")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
sp_i <- read.csv("data/species_3km.csv") %>% filter(Name==sp)
vs.dir <- paste0("vs/", sp_i$Num)
env.in <- readRDS(here(vs.dir, "env_in.rds"))
p <- readRDS(here(vs.dir, "p.rds"))
S <- readRDS(here(vs.dir, "S.rds"))
U <- readRDS(here(vs.dir, "U.rds"))
lam.df <- readRDS(here(vs.dir, "lam_df.rds")) # env & true pop vals
n.cell <- nrow(env.in)



########
## Set sampling details
########
n_samp <- 50  # number of unique samples to average across
O_n <- list(Corr=50, Mech=25)  # number of cells in sample
max_indiv <- list(CA=1e3, IPM=1e2)  # max number of individuals sampled per cell 
O_yr_tmax <- list(Mx=p$tmax, CA=(-2:0)+p$tmax, IPM=p$tmax)  # years to sample
P.i_tmax <- which(lam.df$Surv.S > 5)  # presences: survival past recruit stage
P.pr_tmax <- rep(1, length(P.i_tmax))  # pr(sample cell | presence)
O_yr_tnonEq <- list(Mx=p$tnonEq, CA=(-2:0)+p$tnonEq, IPM=p$tnonEq)
P.i_tnonEq <- which(lam.df$Surv.S_nonEq > 5)
P.pr_tnonEq <- rep(1, length(P.i_tnonEq))
noise <- list(Mx=0.03, # proportion of observed presences that are false
              CA=list(N=0.02,  # N.obs = rnorm(N.true, N.true*N)
                      mu=0.05),  # fec.obs = rnorm(fec.true, fec.true*fec)
              IPM=list(g=0.1,  # sizeNext.obs = rnorm(SizeNext.true, g) 
                       seed=0.05)  # seed.obs = rnorm(seed.true, seed.true*seed)
)


for(s.i in samp.issues) {
  
  ########
  ## Draw samples
  ########
  if(s.i=="nonEq") {
    P.i <- P.i_tnonEq
    P.pr <- P.pr_tnonEq
    O_yr <- O_yr_tnonEq
  } else {
    P.i <- P.i_tmax
    P.pr <- P.pr_tmax
    O_yr <- O_yr_tmax
  }
  if(s.i=="sampBias") P.pr <- env.in$prSamp[P.i]
  
  Corr.sample <- map(1:n_samp, ~sample(P.i, O_n$Corr, replace=F, prob=P.pr))
  Mech.sample <- map(1:n_samp, ~sample(P.i, O_n$Mech, replace=F, prob=P.pr))
  
  O_Mx <- map(Corr.sample, ~(1:n.cell %in% .))
  O_CA <- sample_for_CA(sp, S, lam.df, Mech.sample, O_yr, max_indiv)
  O_IPM <- sample_for_IPM(p, S, lam.df, Mech.sample, O_yr, max_indiv)
  
  
  ########
  ## Impose sampling error
  ########
  if(s.i=="noise") {
    O_Mx <- add_noise_Mx(noise$Mx, Corr.sample, n.cell, P.i)
    O_CA <- add_noise_CA(O_CA, noise$CA)
    O_IPM <- add_noise_IPM(O_IPM, U$lo, U$hi, noise$IPM)
  }
  
  
  ########
  ## Store samples
  ########
  if(overwrite) {
    saveRDS(O_Mx, here(vs.dir, paste0("O_Mx_", s.i, ".rds")))
    saveRDS(O_CA, here(vs.dir, paste0("O_CA_", s.i, ".rds")))
    saveRDS(O_IPM, here(vs.dir, paste0("O_IPM_", s.i, ".rds")))
  }
}




