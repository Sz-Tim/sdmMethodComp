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
sampling.issue <- c("none", "noise", "geog", "bias")[4]

# load workspace
pkgs <- c("tidyverse", "magrittr")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
source("code/fn_IPM.R"); source("code/fn_aux.R"); source("code/fn_sim.R")
env.in <- readRDS(paste0("out/", sp, "_env_in.rds"))
n.cell <- nrow(env.in)
p <- readRDS(paste0("out/", sp, "_p.rds"))
S <- readRDS(paste0("out/", sp, "_S.rds"))
U <- readRDS(paste0("out/", sp, "_U.rds"))
lam.df <- readRDS(paste0("out/", sp, "_lam_df.rds")) # env + true pop values


########
## Set sampling details
########
n_samp <- 5  # number of unique samples to average across
O_n <- list(Corr=100, Mech=20)  # number of cells in sample
O_yr <- list(Mx=p$tmax, CA=(p$tmax-10):p$tmax, IPM=p$tmax)  # years to sample
P.i <- which(lam.df$Surv.S > 0)  # presences: survival past recruit stage
P.pr <- rep(1, length(P.i))  # pr(sample cell | presence)
prop.sampled <- 1  # proportion of individuals sampled per sampled cell 
geog.excl <- which(env.in$x > 20 & env.in$y > 45)
noise <- list(Mx=0.2, # proportion of observed presences that are false
            CA=NA,
            IPM=list(s=0.1,  # proportion of incorrectly assessed surv
                     g=0.2,  # sizeNext.obs = SizeNext.true + rnorm(0,SD) 
                     fl=0.1,  # proportion of incorrectly assessed fl
                     seed=10)  # seed.obs = seed.true + rnorm(0,SD)
            )


########
## Draw samples
########
# determine cells
if(sampling.issue=="geog") {
  P.pr[P.i %in% geog.excl] <- 0
} else if(sampling.issue=="bias") {
  P.pr <- P.pr * env.in$rdLen[P.i]
}
Corr.sample <- map(1:n_samp, ~sample(P.i, O_n$Corr, replace=F, prob=P.pr))
Mech.sample <- map(1:n_samp, ~sample(P.i, O_n$Mech, replace=F, prob=P.pr))

# MaxEnt: presences
O_Mx <- map(Corr.sample, ~(1:nrow(lam.df) %in% .))

# CA: population samples
O_CA <- vector("list", n_samp)
for(s in 1:n_samp) {
  CA.d <- vector("list", O_n$Mech)
  for(i in Mech.sample[[s]]) {
    CA.d[[i]] <- data.frame(S$d[[i]]) %>% filter(yr %in% O_yr$CA) %>%
      add_column(id.inbd=i) %>%
      full_join(env.in[i,], by="id.inbd")
    CA.d[[i]] <- sample_frac(CA.d[[i]], prop.sampled)
  }
  O_CA[[s]] <- do.call(rbind, CA.d)
}

# IPM: population samples
O_IPM <- vector("list", n_samp)
for(s in 1:n_samp) {
  IPM.d <- vector("list", O_n$Mech)
  for(i in Mech.sample[[s]]) {
    IPM.d[[i]] <- data.frame(S$d[[i]]) %>% filter(yr %in% O_yr$IPM) %>%
      mutate(size2=size^2, size3=size^3) %>%
      add_column(id.inbd=i) %>%
      full_join(env.in[i,], by="id.inbd")
    IPM.d[[i]] <- sample_frac(IPM.d[[i]], prop.sampled)
  }
  O_IPM[[s]] <- do.call(rbind, IPM.d)
}

# add error
if(sampling.issue=="noise") {
  # MaxEnt: substitute in some % false presences
  n.err <- noise$Mx*O_n$Corr
  for(s in 1:n_samp) {
    Corr.sample[[s]][1:n.err] <- sample((1:n.cell)[-P.i], n.err, F)
  }
  O_Mx <- map(Corr.sample, ~(1:nrow(lam.df) %in% .))
  # CA: add noise... 
  # IPM: add noise... 
  for(s in 1:n_samp) {
    n_obs <- nrow(O_IPM[[s]])
    # adjust survival
    # adjust growth, flowering based on survival
    O_IPM[[s]]$sizeNext %<>% pmin(. + rnorm(n_obs, 0, noise$IPM$g), U$hi)
    O_IPM[[s]]$sizeNext %<>% pmax(., U$lo)
    # adjust flowering
    # adjust seed based on flowering
    O_IPM[[s]]$seed %<>% pmax(. + rnorm(n_obs, 0, noise$IPM$seed), 0)
  }
}


########
## Store samples
########
saveRDS(O_Mx, paste0("out/", sp, "_O_Mx_", sampling.issue, ".rds"))
saveRDS(O_CA, paste0("out/", sp, "_O_CA_", sampling.issue, ".rds"))
saveRDS(O_IPM, paste0("out/", sp, "_O_IPM_", sampling.issue, ".rds"))



