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
spp.virt <- c(barberry="sp1", lindera="sp2", 
              garlic_mustard="sp3", tower_mustard="sp4")
sp <- names(spp.virt)[1]
overwrite <- TRUE
env.f <- "data/ENF_10km.csv"  # file with environmental data
clim_X <- paste0("bio10_", c(1, 12))
Sys.setenv("MC_CORES"=4)

# load workspace
pkgs <- c("gbPopMod", "tidyverse", "magrittr", "here", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
agg.sp <- read.csv(here("data/PNAS_2017/", ifelse(grepl("mustard", sp),
                                                  "aggLC_mustard.csv", 
                                                  "aggLC_woody.csv")))
L <- build_landscape(f=here(env.f), nlcd_agg=agg.sp, clim_X=clim_X,
                     x_max=Inf, y_max=150) 
n.cell <- sum(L$env.rct$inbd)



########
## Set species traits
########
# p=list(n=30,  # ncells in IPM matrix
#        tmax=100,  # time steps for NDD & simulations
#        n0=100,  # initial pop sizes
#        prop_init=0.01,  # proportion of cells with initial populations
#        z.rng=c(1,12),  # initial size range
#        s_z=c(-4, 2.1, -.09),  # b1 + b2*z + b3*z^2
#        s_x=c(-1, -0.3, -2, 0, -1, -0.2),  # b1*x1 + ...
#        g_z=c(.2, 2, -0.1),  # b1 + b2*z + b3*z^2
#        g_x=c(-0.2, 0, -2, 0.2, -1, 0.3, -1, 0, -3, 1, 0.5, 3),  # b1*x1 + ...
#        g_sig=1,  # growth ~ N(E, g_sig)
#        fl_z=c(-1.5, .1, .1),  # b1 + b2*z + b3*z^2
#        fl_x=c(-1, 0, 0, -1, -1, 0.2, -1, -0.1),  # b1*x1 + ...
#        seed_z=c(2, 0.5, -.03),  # b1 + b2*z + b2*z^2
#        seed_x=c(-0.7, -0.2, -1, -0.1, -2),  # b1*x1 + ...
#        rcr_z=c(1.5, 0.4),  # N(mean=rcrt1, sd=rcrt2)
#        p_est=0.03,  # p(establishment)
#        NDD=T,  # negative density dependent p_est
#        rcr_SB=0.5,  # p(recruit from seedbank)
#        rcr_dir=0,  # p(recruit directly)
#        s_SB=0.75,  # p(survive in seedbank additional year)
#        sdd_max=4,  # max SDD distance in cells
#        sdd_rate=3,  # SDD dispersal rate
#        ldd=1,  # number of LDD events per year
#        bird_hab=c(1,1,1,1,1)  # bird habitat preferences among LC types
# )
p <- fit_PNAS_species(sp=sp, nlcd_agg=agg.sp, clim_X=clim_X, n_z=3)
p$n <- 30
p$tmax <- 50
p$n0 <- 100
p$prop_init <- 0.01
p$NDD <- T
p$sdd_max <- 4
p$sdd_rate <- 3
p$ldd <- 1
p$bird_hab <- c(1,1,1,1,1)
p$NDD_n <- p$n0/10  # mean number of recruits if NDD
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
sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=L$env.rct.unscaled, 
                        lc.col=tail(1:ncol(L$env.rct.unscaled), 
                                    n_distinct(agg.sp$agg)),
                        g.p=list(sdd.max=p$sdd_max, 
                                 sdd.rate=p$sdd_rate, 
                                 bird.hab=p$bird_hab))
sdd.j <- mclapply(1:n.cell, function(x) which(sdd.pr$i[,,2,]==L$env.in$id[x], 
                                              arr.ind=T))
p.ij <- mclapply(1:n.cell, function(x) sdd.pr$i[,,1,][sdd.j[[x]]]) 
# NOTE: sdd.pr[,,2,] indexes based on `id` (id for each cell in grid) instead  
# of `id.inbd` (id for inbound cells only), but sdd.pr[,,,i] includes only
# inbound cells, so the layer index aligns with `id.inbd`. This makes 
# identifying much simpler, but requires looking up the corresponding id's


########
## Generate underlying IPM and simulated data
########
# Initial populations
N_init <- rep(0, n.cell)
# N_init[sample.int(n.cell, p$prop_init*n.cell, replace=F)] <- p$n0
N_init[sample(filter(L$env.in, x>35 & y<20)$id.inbd, 
              p$prop_init*n.cell, replace=F)] <- p$n0

# Use assigned slopes to fill IPM matrix
U <- fill_IPM_matrices(n.cell, buffer=0.75, discrete=1, p, n_z, n_x, 
                       X, sdd.j, p.ij, N_init, verbose=T)

# Ground Truth: generate simulated data
S <- simulate_data(n.cell, U$lo, U$hi, p, X, n_z, sdd.pr$i, sdd.j, N_init, T)
# Aggregate results
lam.df <- L$env.in %>%
  mutate(lambda=apply(U$IPMs, 3, function(x) Re(eigen(x)$values[1])),
         lam.S=map_dbl(S$d, ~sum(.$surv[.$yr==p$tmax], na.rm=T))/
           (map_dbl(S$d, ~sum(.$surv[.$yr==(p$tmax-1)], na.rm=T))+.01),
         nSeed=S$nSd[,p$tmax], 
         D=S$D[,p$tmax], 
         B0=S$B[,1], 
         Btmax=S$B[,p$tmax+1], 
         N.S=map_dbl(S$d, ~sum(!is.na(.$sizeNext[.$yr==p$tmax]))),
         Surv.S=map_dbl(S$d, ~sum(.$surv[.$yr==p$tmax], na.rm=T)),
         Rcr.S=map_dbl(S$d, ~sum(is.na(.$size[.$yr==p$tmax]))),
         nSdStay=nSeed*(1-p$p_emig), 
         nSdLeave=nSeed*p$p_emig,
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
  if(!dir.exists(here("vs", spp.virt[[sp]]))) {
    dir.create(here("vs", spp.virt[[sp]]), recursive=T)
  }
  saveRDS(L$scale.i, here("vs", spp.virt[[sp]], "cov_scale.rds"))
  saveRDS(L$env.rct, here("vs", spp.virt[[sp]], "env_rct.rds"))
  saveRDS(L$env.rct.unscaled, here("vs", spp.virt[[sp]], "env_rct_unscaled.rds"))
  saveRDS(L$env.in, here("vs", spp.virt[[sp]], "env_in.rds"))
  saveRDS(p, here("vs", spp.virt[[sp]], "p.rds"))
  saveRDS(N_init, here("vs", spp.virt[[sp]], "N_init.rds"))
  saveRDS(sdd.pr, here("vs", spp.virt[[sp]], "sdd.rds"))
  saveRDS(sdd.j, here("vs", spp.virt[[sp]], "sdd_j.rds"))
  saveRDS(p.ij, here("vs", spp.virt[[sp]], "p_ij.rds"))
  saveRDS(U, here("vs", spp.virt[[sp]], "U.rds"))
  saveRDS(S, here("vs", spp.virt[[sp]], "S.rds"))
  saveRDS(lam.df, here("vs", spp.virt[[sp]], "lam_df.rds"))
}



