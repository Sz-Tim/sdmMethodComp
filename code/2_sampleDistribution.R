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
sampling.issue <- c("none", "noise", "geogBias", "sampBias")[2]

# load workspace
pkgs <- c("tidyverse", "magrittr", "here")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
env.in <- readRDS(here(paste0("out/", sp, "_env_in.rds")))
n.cell <- nrow(env.in)
p <- readRDS(here(paste0("out/", sp, "_p.rds")))
S <- readRDS(here(paste0("out/", sp, "_S.rds")))
U <- readRDS(here(paste0("out/", sp, "_U.rds")))
lam.df <- readRDS(here(paste0("out/", sp, "_lam_df.rds"))) # env + true pop vals


########
## Set sampling details
########
n_samp <- 10  # number of unique samples to average across
O_n <- list(Corr=100, Mech=25)  # number of cells in sample
O_yr <- list(Mx=p$tmax, CA=(-1:0)+p$tmax, IPM=p$tmax)  # years to sample
P.i <- which(lam.df$Surv.S > 5)  # presences: survival past recruit stage
P.pr <- rep(1, length(P.i))  # pr(sample cell | presence)
prop.sampled <- 1  # proportion of individuals sampled per sampled cell 
geog.excl <- which(env.in$y < 30)
noise <- list(Mx=0.2, # proportion of observed presences that are false
            CA=list(N=0.02,  # N.obs = rnorm(N.true, N.true*N)
                    fec=0.05),  # fec.obs = rnorm(fec.true, fec.true*fec)
            IPM=list(s=0,  # proportion of incorrectly assessed surv
                     g=0.1,  # sizeNext.obs = rnorm(SizeNext.true, g) 
                     fl=0,  # proportion of incorrectly assessed fl
                     seed=0.05)  # seed.obs = rnorm(seed.true, seed.true*seed)
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
set.seed(11)
Corr.sample <- map(1:n_samp, ~sample(P.i, O_n$Corr, replace=F, prob=P.pr))
Mech.sample <- map(1:n_samp, ~sample(P.i, O_n$Mech, replace=F, prob=P.pr))

# MaxEnt: presences
O_Mx <- map(Corr.sample, ~(1:nrow(lam.df) %in% .))

# CA: population samples
O_CA <- vector("list", n_samp)
for(s in 1:n_samp) {
  CA.d <- CA.B <- CA.D <- vector("list", O_n$Mech)
  for(j in seq_along(Mech.sample[[s]])) {
    i <- Mech.sample[[s]][j]
    CA.d[[j]] <- data.frame(S$d[[i]]) %>% filter(yr %in% O_yr$CA)
    CA.d[[j]] <- sample_frac(CA.d[[j]], prop.sampled)
    CA.d[[j]] <- CA.d[[j]] %>% group_by(yr) %>% 
      summarise(N=sum(!is.na(size)), 
                s.ad.0=sum(surv[age>2]==0, na.rm=TRUE),
                s.ad.1=sum(surv[age>2]==1, na.rm=TRUE),
                s.jv.0=sum(surv[age<3]==0, na.rm=TRUE),
                s.jv.1=sum(surv[age<3]==1, na.rm=TRUE),
                f.0=sum(fl==0, na.rm=TRUE),
                f.1=sum(fl==1, na.rm=TRUE),
                fec=mean(seed, na.rm=TRUE) %>% round,
                nSeed=sum(seed, na.rm=TRUE),
                N.rcr=sum(is.na(size))) %>%
      mutate(fec=ifelse(is.nan(fec), 0, fec)) %>%
      ungroup() %>%
      mutate(lambda=N/lag(N,1)) %>%
      add_column(id.inbd=i) %>%
      full_join(env.in[i,], by="id.inbd")
    CA.B[[j]] <- S$B[i,O_yr$CA]
    CA.D[[j]] <- S$D[i,O_yr$CA]
  }
  O_CA[[s]] <- list(d=do.call(rbind, CA.d),
                    B=do.call(rbind, CA.B),
                    D=do.call(rbind, CA.D))
  O_CA[[s]]$d$p.est <- c(t(with(O_CA[[s]], d$N.rcr/(d$nSeed*p$p_emig*p$rcr_dir + 
                                                  B*p$rcr_SB + D*p$rcr_dir))))
}

# IPM: population samples
O_IPM <- vector("list", n_samp)
for(s in 1:n_samp) {
  IPM.d <- vector("list", O_n$Mech)
  for(j in seq_along(Mech.sample[[s]])) {
    i <- Mech.sample[[s]][j]
    IPM.d[[j]] <- data.frame(S$d[[i]]) %>% filter(yr %in% O_yr$IPM) %>%
      mutate(size2=size^2, size3=size^3) %>%
      add_column(id.inbd=i) %>%
      full_join(env.in[i,], by="id.inbd")
    IPM.d[[j]] <- sample_frac(IPM.d[[j]], prop.sampled)
  }
  O_IPM[[s]] <- do.call(rbind, IPM.d)
}

# add error
if(sampling.issue=="noise") {
  # MaxEnt: substitute false presences
  n.err <- noise$Mx*O_n$Corr
  for(s in 1:n_samp) {
    Corr.sample[[s]][1:n.err] <- sample((1:n.cell)[-P.i], n.err, F)
  }
  O_Mx <- map(Corr.sample, ~(1:nrow(lam.df) %in% .))
  # CA: add error to N and seed counts
  for(s in 1:n_samp) {
    n_obs <- nrow(O_CA[[s]]$d)
    O_CA[[s]]$d <- O_CA[[s]]$d %>% 
      mutate(N=pmax(round(N + rnorm(n_obs, 0, N*noise$CA$N)), 0),
             fec=pmax(round(fec + rnorm(n_obs, 0, fec*noise$CA$fec))), 0) %>%
      group_by(id) %>%
      mutate(lambda=N/lag(N,1))
  }
  # IPM: add error to growth measurements and seed counts
  for(s in 1:n_samp) {
    n_obs <- nrow(O_IPM[[s]])
    O_IPM[[s]]$sizeNext <- O_IPM[[s]]$sizeNext + rnorm(n_obs, 0, noise$IPM$g)
    O_IPM[[s]]$sizeNext <- pmin(O_IPM[[s]]$sizeNext, U$hi)
    O_IPM[[s]]$sizeNext <- pmax(O_IPM[[s]]$sizeNext, U$lo)
    O_IPM[[s]]$seed <- pmax(round(O_IPM[[s]]$seed +
                              rnorm(n_obs, 0, O_IPM[[s]]$seed*noise$IPM$seed)), 
                            0)
  }
}


########
## Store samples
########
if(overwrite) {
  saveRDS(O_Mx, here(paste0("out/", sp, "_O_Mx_", sampling.issue, ".rds")))
  saveRDS(O_CA, here(paste0("out/", sp, "_O_CA_", sampling.issue, ".rds")))
  saveRDS(O_IPM, here(paste0("out/", sp, "_O_IPM_", sampling.issue, ".rds")))
}




