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
env.f <- "data/ENF_3km.csv"  # file with environmental data
clim_X <- paste0("bio10_", c(5, "prMay"))
habitat <- 3
max_z_pow <- 1
n.cores <- 8
x_min <- 675
x_max <- Inf
y_min <- 0
y_max <- 250

# load workspace
pkgs <- c("gbPopMod", "tidyverse", "magrittr", "here", "doSNOW", "foreach")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
nlcd.sp <- read.csv(here("data/PNAS_2017/", ifelse(grepl("mustard", sp),
                                                   "aggLC_mustard.csv", 
                                                   "aggLC_woody.csv")))
L <- build_landscape(f=here(env.f), nlcd_agg=nlcd.sp, clim_X=clim_X,
                     x_min, x_max, y_min, y_max) 
n.cell <- sum(L$env.rct$inbd)



########
## Set species traits
########
p <- fit_PNAS_species(sp, env.f, nlcd.sp, clim_X, FALSE, max_z_pow, habitat,
                      x_min, x_max, y_min, y_max)
p$n <- 50
p$tmax <- 150
p$n0 <- 10
p$prop_init <- 0.001
p$NDD <- T
p$sdd_max <- 7
p$sdd_rate <- 1.4
p$ldd <- 2
p$bird_hab <- c(.32, .36, .05, .09, .09)
p$NDD_n <- p$n0/10  # mean number of recruits if NDD
p$p_emig <- pexp(0.5, p$sdd_rate, lower.tail=F) # p(seed emigrants)
n_z <- list(s=length(p$s_z),  # n size covariates for each vital rate
            g=length(p$g_z),
            fl=length(p$fl_z), 
            seed=length(p$seed_z))
n_x <- list(s=length(p$s_x), # n env covariates for each vital rate
            g=length(p$g_x),
            fl=length(p$fl_x), 
            seed=length(p$seed_x),
            germ=length(p$germ_x))
X <- map(n_x, ~as.matrix(L$env.in[,1:.]))  # env covariates for each vital rate
if(!is.null(X$germ)) X$germ <- cbind(1, X$germ[,-n_x$germ])
# p$p_emig <- 0
sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=L$env.rct.unscaled,
                        lc.col=tail(1:ncol(L$env.rct.unscaled),
                                    n_distinct(nlcd.sp$agg)),
                        g.p=list(sdd.max=p$sdd_max,
                                 sdd.rate=p$sdd_rate,
                                 bird.hab=p$bird_hab))
sdd.df <- data.frame(i=rep(1:n.cell, times=map_int(sdd.pr$sp, length)),
                     j=unlist(map(sdd.pr$sp, ~as.numeric(names(.)))),
                     pr=unlist(sdd.pr$sp))
sdd.j <- map(L$env.in$id, ~sdd.df$i[sdd.df$j==.])
p.ij <- map(L$env.in$id, ~sdd.df$pr[sdd.df$j==.])
sdd.df$j_in <- unlist(sdd.j)
# NOTE: sdd.pr[,,2,] indexes based on `id` (id for each cell in grid) instead  
# of `id.inbd` (id for inbound cells only), but sdd.pr[,,,i] includes only
# inbound cells, so the layer index aligns with `id.inbd`. This makes 
# identifying much simpler, but requires looking up the corresponding id's


########
## Generate underlying IPM and simulated data
########
# Initial populations
N_init <- rep(0, n.cell)
N_init[sample(filter(L$env.in, x>725 & x<765 & y>200)$id.inbd,
              p$prop_init*n.cell, replace=F)] <- p$n0

# Use assigned slopes to fill IPM matrix
U <- fill_IPM_matrices(n.cell, buffer=0.1, discrete=1, p, n_z, n_x, 
                       X, sdd.j, p.ij, verbose=T)

# Ground Truth: generate simulated data
S <- simulate_data(n.cell, U$lo, U$hi, p, X, n_z, sdd.df, p.ij, N_init, 
                   save_yrs=(-2:0)+p$tmax, T)

# Aggregate results
lam.df <- L$env.in %>%
  mutate(lambda=apply(U$IPMs, 3, function(x) Re(eigen(x)$values[1])),
         lam.S=map_dbl(S$d, ~sum(.$surv[.$yr==p$tmax], na.rm=T))/
           (map_dbl(S$d, ~sum(.$surv[.$yr==(p$tmax-1)], na.rm=T))+.01),
         nSeed=S$nSd[,dim(S$nSd)[2]], 
         D=S$D[,dim(S$D)[2]], 
         B=S$B[,dim(S$B)[2]], 
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
         max.age=map_dbl(S$d, ~max(.$age)),
         pr_Immigrant=map_dbl(p.ij, sum),
         s=antilogit(c(cbind(1, mean(p$z.rng)/2, X$s) %*% c(p$s_z, p$s_x))),
         g=c(cbind(1, mean(p$z.rng), X$g) %*% c(p$g_z, p$g_x)),
         germ=antilogit(c(X$germ %*% p$germ_x)))

library(viridis)
lam.gg <- ggplot(lam.df, aes(x=lon, y=lat)) + 
  scale_fill_viridis(name="", option="B") +
  ggtitle(paste0(sp, ": 3km x 3km, favorable habitat"))
lam.gg + geom_tile(aes(fill=log(lambda))) + labs(subtitle="log(lambda)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=log(Surv.S))) + labs(subtitle="log(N)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=log(nSeed))) + labs(subtitle="log(Seed production)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=nSeed/Surv.S)) + labs(subtitle="Per capita seed production") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=log(round(D)))) + labs(subtitle="log(Immigrant seeds)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=log(nSdStay+round(D)))) + labs(subtitle="log(Propagule pressure)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=log(B))) + labs(subtitle="log(Seed bank)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=lambda>1)) +  labs(subtitle="lambda > 1") +
  scale_fill_manual("", values=c("gray30", "dodgerblue")) +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=Surv.S>0)) + labs(subtitle="N > 0") +
  scale_fill_manual("", values=c("gray30", "dodgerblue")) +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=s)) + labs(subtitle="s: mean(z.rng)/2") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=g)) + labs(subtitle="g = mn(growth): mean(z.rng)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
lam.gg + geom_tile(aes(fill=germ)) + labs(subtitle="pr(germination)") +
  geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  



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



