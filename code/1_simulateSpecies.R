# 1: Simulate species
# Comparison of SDM approaches
# Tim Szewczyk

# This script generates data for hypothetical species. The simulation process is
# based on an IPM framework. Both the underlying IPM objects and the simulated
# individuals are produced and stored. This serves as the perfectly known true
# distribution from which samples are taken for fitting each SDM.

########
## Setup
########
# file specifications
sp <- "sp1"
overwrite <- TRUE
env.f <- "data/landcover_5km.csv"  # file with environmental data

# load workspace
pkgs <- c("gbPopMod", "tidyverse", "magrittr", "here")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))

L <- build_landscape(f=here(env.f), 
                     x_max=Inf, # ncol in landscape; Inf for full dataset
                     y_max=Inf) # nrow in landscape; Inf for full dataset
n.cell <- sum(L$env.rct$inbd)


########
## Set species traits
########
p=list(n=30,  # ncells in IPM matrix
       tmax=30,  # time steps for NDD & simulations
       n0=100,  # initial pop sizes
       prop_init=0.2,  # proportion of cells with initial populations
       z.rng=c(1,12),  # initial size range
       s_z=c(-8, 2.1, -.09),  # b1 + b2*z + b3*z^2
       s_x=c(3, -.1, -2, -.1, 2, -2, -.4),  # b1*x1 + ...
       g_z=c(.2, 2, -0.1),  # b1 + b2*z + b3*z^2
       g_x=c(2, -.1, 2, -.1, 2, -2),  # b1*x1 + ...
       g_sig=1,  # growth ~ N(E, g_sig)
       fl_z=c(-1.5, .1, .1),  # b1 + b2*z + b3*z^2
       fl_x=c(-2, -.1, -2, -.1, 1, 1),  # b1*x1 + ...
       seed_z=c(2, 0.5, -.03),  # b1 + b2*z + b2*z^2
       seed_x=c(1, -.1, -1, -.1, .2),  # b1*x1 + ...
       rcr_z=c(1.5, 0.4),  # N(mean=rcrt1, sd=rcrt2)
       p_est=0.03,  # p(establishment)
       NDD=T,  # negative density dependent p_est
       rcr_SB=0.8,  # p(recruit from seedbank)
       rcr_dir=0.1,  # p(recruit directly)
       s_SB=0.8,  # p(survive in seedbank additional year)
       sdd_max=5,  # max SDD distance in cells
       sdd_rate=1,  # SDD dispersal rate
       bird_hab=c(1,1,1,1,1)  # bird habitat preferences among LC types
)
p$NDD_n <- p$n0/3  # mean number of recruits if NDD
p$p_emig <- pexp(0.5, p$sdd_rate, lower.tail=F) # p(seed emigrants)
n_z <- list(s=length(p$s_z),  # n size covariates for each vital rate
            g=length(p$g_z),
            fl=length(p$fl_z), 
            seed=length(p$seed_z))
n_x <- list(s=length(p$s_x), # n env covariates for each vital rate
            g=length(p$g_x),
            fl=length(p$fl_x), 
            seed=length(p$seed_x))
X <- map(n_x, ~as.matrix(L$env.in[,1:.]))  # env covariates for each vital rate
sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=L$env.rct.unscaled, lc.col=8:12,
                        g.p=list(sdd.max=p$sdd_max, 
                                 sdd.rate=p$sdd_rate, 
                                 bird.hab=p$bird_hab))
# NOTE: sdd.pr[,,2,] indexes based on `id` (id for each cell in grid) instead  
# of `id.inbd` (id for inbound cells only), but sdd.pr[,,,i] includes only
# inbound cells, so the layer index aligns with `id.inbd`. This makes 
# identifying much simpler, but requires looking up the corresponding id's


########
## Generate underlying IPM and simulated data
########
# Initial populations
N_init <- rep(0, n.cell)
N_init[sample.int(n.cell, p$prop_init*n.cell, replace=F)] <- p$n0
N_init[arrange(lam.df, desc(Surv.S))$id.inbd[1:(p$prop_init*n.cell)]] <- p$n0

# Use assigned slopes to fill IPM matrix
U <- fill_IPM_matrices(n.cell, buffer=0.75, discrete=1, p, n_z, n_x, 
                       X, sdd.pr, L$env.in$id, N_init)

# Ground Truth: generate simulated data
S <- simulate_data(n.cell, U$lo, U$hi, p, X, n_z, sdd.pr, U$sdd.j, N_init)

# Aggregate results
lam.df <- L$env.in %>%
  mutate(lambda=apply(U$IPMs, 3, function(x) Re(eigen(x)$values[1])),
         lam.S=S$N_sim[,p$tmax]/(S$N_sim[,p$tmax-1]+.01),
         nSeed=S$nSd[,p$tmax], 
         D=S$D[,p$tmax], 
         B0=S$B[,1], 
         Btmax=S$B[,p$tmax+1], 
         N.S=map_dbl(S$d, ~sum(!is.na(.$sizeNext[.$yr==p$tmax]))),
         Surv.S=map_dbl(S$d, ~sum(.$surv[.$yr==p$tmax], na.rm=T)),
         Rcr.S=map_dbl(S$d, ~sum(is.na(.$size[.$yr==p$tmax]))),
         nSdStay=nSeed*(1-p$p_emig), 
         nSdLeave=nSeed*p$p_emig,
         N.U=apply(U$Nt[-1,,p$tmax+1],2,sum), 
         lam.U=U$lam.t[,p$tmax], 
         mn.age.z=map_dbl(S$d, ~mean(.$age[.$yr==p$tmax & 
                                           !is.na(.$sizeNext) &
                                           !is.na(.$size)])),
         med.age.z=map_dbl(S$d, ~median(.$age[.$yr==p$tmax & 
                                              !is.na(.$sizeNext) &
                                              !is.na(.$size)])),
         mn.age=map_dbl(S$d, ~mean(.$age[.$yr==p$tmax])),
         med.age=map_dbl(S$d, ~median(.$age[.$yr==p$tmax])),
         max.age=map_dbl(S$d, ~max(.$age)))


########
## Store true species distribution
########
if(overwrite) {
  saveRDS(L$scale.i, here(paste0("out/", sp, "_cov_scale.rds")))
  saveRDS(L$env.rct, here(paste0("out/", sp, "_env_rct.rds")))
  saveRDS(L$env.rct.unscaled, here(paste0("out/", sp, "_env_rct_unscaled.rds")))
  saveRDS(L$env.in, here(paste0("out/", sp, "_env_in.rds")))
  saveRDS(p, here(paste0("out/", sp, "_p.rds")))
  saveRDS(N_init, here(paste0("out/", sp, "_N_init.rds")))
  saveRDS(sdd.pr, here(paste0("out/", sp, "_sdd.rds")))
  saveRDS(U, here(paste0("out/", sp, "_U.rds")))
  saveRDS(S, here(paste0("out/", sp, "_S.rds")))
  saveRDS(lam.df, here(paste0("out/", sp, "_lam_df.rds")))
}



