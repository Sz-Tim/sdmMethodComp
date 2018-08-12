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
sp <- "sp1"
overwrite <- TRUE
sampling.issue <- c("none", "noise", "geogBias", "sampBias")[4]

# load workspace
pkgs <- c("tidyverse", "magrittr", "here")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), ~source(here(.)))
env.in <- readRDS(here("vs", sp, "env_in.rds"))
p <- readRDS(here("vs", sp, "p.rds"))
S <- readRDS(here("vs", sp, "S.rds"))
U <- readRDS(here("vs", sp, "U.rds"))
lam.df <- readRDS(here("vs", sp, "lam_df.rds")) # env & true pop vals
n.cell <- nrow(env.in)


########
## Set sampling details
########
n_samp <- 5  # number of unique samples to average across
O_n <- list(Corr=50, Mech=25)  # number of cells in sample
O_yr <- list(Mx=p$tmax, CA=(-2:0)+p$tmax, IPM=p$tmax)  # years to sample
P.i <- which(lam.df$Surv.S > 5)  # presences: survival past recruit stage
P.pr <- rep(1, length(P.i))  # pr(sample cell | presence)
prop.sampled <- 1  # proportion of individuals sampled per sampled cell 
geog.excl <- which(env.in$x < 20)
noise <- list(Mx=0.2, # proportion of observed presences that are false
              CA=list(N=0.02,  # N.obs = rnorm(N.true, N.true*N)
                      mu=0.05),  # fec.obs = rnorm(fec.true, fec.true*fec)
              IPM=list(g=0.1,  # sizeNext.obs = rnorm(SizeNext.true, g) 
                       seed=0.05)  # seed.obs = rnorm(seed.true, seed.true*seed)
)


########
## Impose sampling bias
########
if(sampling.issue=="geogBias") {
  P.pr[P.i %in% geog.excl] <- 0
} else if(sampling.issue=="sampBias") {
  P.pr <- P.pr * env.in$rdLen[P.i]
}


########
## Draw samples
########
Corr.sample <- map(1:n_samp, ~sample(P.i, O_n$Corr, replace=F, prob=P.pr))
Mech.sample <- map(1:n_samp, ~sample(P.i, O_n$Mech, replace=F, prob=P.pr))

O_Mx <- map(Corr.sample, ~(1:n.cell %in% .))
O_CA <- sample_for_CA(S, Mech.sample, O_n, O_yr, n_samp, prop.sampled)
O_IPM <- sample_for_IPM(S, Mech.sample, O_n, O_yr, n_samp, prop.sampled)


########
## Impose sampling error
########
if(sampling.issue=="noise") {
  O_Mx <- add_noise_Mx(O_Mx, noise$Mx, O_n, n_samp, Corr.sample, n.cell, P.i)
  O_CA <- add_noise_CA(O_CA, noise$CA, n_samp)
  O_IPM <- add_noise_IPM(O_IPM, U$lo, U$hi, noise$IPM, n_samp)
}


########
## Store samples
########
if(overwrite) {
  saveRDS(O_Mx, here("vs", sp, paste0("O_Mx_", sampling.issue, ".rds")))
  saveRDS(O_CA, here("vs", sp, paste0("O_CA_", sampling.issue, ".rds")))
  saveRDS(O_IPM, here("vs", sp, paste0("O_IPM_", sampling.issue, ".rds")))
}




