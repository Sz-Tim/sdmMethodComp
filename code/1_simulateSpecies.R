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
sp <- c("barberry", "garlic_mustard")[1]
res <- "10km" # "3km", "5km", "10km", "50km"
overwrite <- TRUE
plots <- TRUE
clim_X <- paste0("bio10_", c(6, "prMay"))
habitat <- 4
max_z_pow <- 1
n.cores <- 24
x_min <- 0#200#675#
x_max <- Inf
y_min <- 0
y_max <- Inf#75#250#

# load workspace
pkgs <- c("gbPopMod", "tidyverse", "magrittr", "here", "doSNOW", "foreach")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
env.f <- paste0("data/ENF_", res, ".csv")
sp_i <- read.csv(paste0("data/species_", res, ".csv")) %>% filter(Name==sp)
nlcd.sp <- read.csv(here("data/PNAS_2017/", sp_i$LC_f))
L <- build_landscape(env.f, nlcd.sp, x_min, x_max, y_min, y_max) 
n.cell <- sum(L$env.rct$inbd)



########
## Set species traits
########
p <- fit_PNAS_species(sp, env.f, nlcd.sp, clim_X, FALSE, max_z_pow, habitat,
                      x_min, x_max, y_min, y_max)
p$s_x <- c(-2.75, -1.25, 3.5, 0)
p$g_x <- c(-1.25, -0.75, -0.3, -0.3)
p$germ_x <- c(1.5, -4, -1.75, -2, 0)
p$n <- 20
p$tmax <- 100
p$tnonEq <- floor(p$tmax/3)
p$n0 <- 10
p$prop_init <- 0.001
p$NDD <- T
p$sdd_max <- sp_i$sdd_max
p$sdd_rate <- sp_i$sdd_rate
p$ldd <- sp_i$ldd
p$bird_hab <- c(.32, .36, .05, .09, .09)
p$NDD_n <- 20  # mean number of recruits if NDD
p$K_max <- 1e3  # max abundance for CAd
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
X <- map(n_x, ~as.matrix(L$env.in[,grep(paste(clim_X, collapse="|"), 
                                        names(L$env.in))]))  # env covs
if(!is.null(X$germ)) X$germ <- cbind(1, X$germ[,-n_x$germ])
# p$p_emig <- 0
sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=L$env.rct.unscaled,
                        lc.col=tail(1:ncol(L$env.rct.unscaled),
                                    n_distinct(nlcd.sp$agg)),
                        g.p=list(sdd.max=p$sdd_max,
                                 sdd.rate=p$sdd_rate,
                                 bird.hab=p$bird_hab))
# df row indexes for all cells dispersing to each [[j]]
p.c <- makeCluster(n.cores); registerDoSNOW(p.c)
sdd.ji.rows <- foreach(x=1:n.cell) %dopar% { which(sdd.pr$sp.df$j.idin==x) }
stopCluster(p.c)
# id.in indexes & probabilities for all cells dispersing INTO [[j]]
sdd.ji <- lapply(sdd.ji.rows, function(x) sdd.pr$sp.df$i.idin[x]) 
p.ji <- lapply(sdd.ji.rows, function(x) sdd.pr$sp.df$pr[x]) 
# NOTE: sdd.pr[,,2,] indexes based on `id` (id for each cell in grid) instead  
# of `id.in` (id for inbound cells only), but sdd.pr[,,,i] includes only
# inbound cells, so the layer index aligns with `id.in`. This makes 
# identifying much simpler, but requires looking up the corresponding id's


########
## Generate underlying IPM and simulated data
########
# Initial populations
N_init <- rep(0, n.cell)
# 725 765 200
# N_init[sample(filter(L$env.in, x>200 & y<100)$id.in,
N_init[sample(filter(L$env.in, x>215 & x<230 & y>50 & y<75)$id.in,
              p$prop_init*n.cell, replace=F)] <- p$n0

# Use assigned slopes to fill IPM matrix
U <- fill_IPM_matrices(n.cell, buffer=0, discrete=1, p, n_z, n_x, 
                       X, sdd.ji, p.ji, sp, verbose=T)
if(sp=="garlic_mustard") {
  U$lambda <- sapply(1:n.cell, function(x) iter_lambda(p, U$Ps[,,x], U$Fs[,,x]))
} else {
  U$lambda <- apply(U$IPMs, 3, function(x) Re(eigen(x)$values[1]))
}

# Ground Truth: generate simulated data
S <- simulate_data(n.cell, U$lo, U$hi, p, X, n_z, sdd.ji, p.ji, N_init, sp, 
                   save_yrs=c(p$tnonEq+(-5:0), p$tmax+(-5:0)), T)

# Aggregate results
lam.df <- L$env.in %>%
  mutate(lambda=U$lambda,
         nSeed=S$nSd[,dim(S$nSd)[2]], 
         D=S$D[,dim(S$D)[2]], 
         B=S$B[,dim(S$B)[2]], 
         Surv.S=map_dbl(S$d, ~sum(.$surv[.$yr==p$tmax], na.rm=T)),
         Surv.S_nonEq=map_dbl(S$d, ~sum(.$surv[.$yr==p$tnonEq], na.rm=T)),
         nRepro=map_dbl(S$d, ~sum(.$fl[.$yr==p$tmax], na.rm=T)),
         Rcr.S=map_dbl(S$d, ~sum(is.na(.$size[.$yr==p$tmax]))),
         nSdStay=nSeed*(1-p$p_emig), 
         nSdLeave=nSeed*p$p_emig,
         min.age.fl=map_dbl(S$d, ~min(.$age[.$yr==p$tmax &
                                             .$fl==1]), na.rm=T),
         mn.age.z=map_dbl(S$d, ~mean(.$age[.$yr==p$tmax & 
                                           !is.na(.$sizeNext) &
                                           !is.na(.$size)])),
         med.age.z=map_dbl(S$d, ~median(.$age[.$yr==p$tmax & 
                                              !is.na(.$sizeNext) &
                                              !is.na(.$size)])),
         mn.age=map_dbl(S$d, ~mean(.$age[.$yr==p$tmax])),
         med.age=map_dbl(S$d, ~median(.$age[.$yr==p$tmax])),
         max.age=map_dbl(S$d, ~max(.$age)),
         pr_Immigrant=map_dbl(p.ji, sum),
         s=antilogit(c(cbind(1, mean(p$z.rng)/2, X$s) %*% c(p$s_z, p$s_x))),
         s.max=antilogit(c(cbind(1, max(p$z.rng), X$s) %*% c(p$s_z, p$s_x))),
         g=c(cbind(1, mean(p$z.rng), X$g) %*% c(p$g_z, p$g_x)),
         germ=antilogit(c(X$germ %*% p$germ_x)))

library(viridis)
lam.gg <- ggplot(lam.df, aes(x=lon, y=lat)) + theme_bw() +
  theme(axis.text=element_blank()) + labs(x="", y="") +
  scale_fill_viridis(name="", option="B") +
  ggtitle(paste0(sp, ": 3km x 3km, favorable habitat"))

if(plots) {
  lam.gg + geom_tile(aes(fill=log(lambda))) + 
    labs(subtitle="log(lambda)") +
    geom_point(data=lam.df[which(N_init>0),], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=log(Surv.S))) + 
    labs(subtitle=paste("log(N): year", p$tmax)) +
    geom_point(data=lam.df[which(N_init>0),], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=log(Surv.S_nonEq))) + 
    labs(subtitle=paste("log(N): year", p$tnonEq)) +
    geom_point(data=lam.df[which(N_init>0),], colour="white", shape=1)
  
  lam.gg + geom_tile(aes(fill=log(nSeed))) + 
    labs(subtitle="log(Seed production)") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=nSeed/nRepro)) + 
    labs(subtitle="Per capita seed production") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=log(round(D)))) + 
    labs(subtitle="log(Immigrant seeds)") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=log(nSdStay+round(D)))) + 
    labs(subtitle="log(Propagule pressure)") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=log(B))) + 
    labs(subtitle="log(Seed bank)") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  
  ggplot() + theme_bw() +
    geom_tile(data=lam.df, aes(x=lon, y=lat), fill="gray90") +
    geom_tile(data=filter(lam.df, lambda>1), aes(x=lon, y=lat), 
              fill="darkblue", alpha=0.5) +
    geom_tile(data=filter(lam.df, Surv.S>0), aes(x=lon, y=lat), 
              fill="red", alpha=0.5) +
    theme(axis.text=element_blank()) + labs(x="", y="") +
    ggtitle(paste0(sp, ": 3km x 3km, favorable habitat"))
  
  lam.gg + geom_tile(aes(fill=s)) + 
    labs(subtitle="s: mean(z.rng)/2") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=s.max)) + 
    labs(subtitle="s: max(z.rng)") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=g)) + 
    labs(subtitle="g = mn(growth): mean(z.rng)") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=germ)) + 
    labs(subtitle="pr(germination)") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
}




########
## Store true species distribution
########
if(overwrite) {
  vs.dir <- paste0("vs/", sp_i$Num)
  if(!dir.exists(here(vs.dir))) dir.create(here(vs.dir), recursive=T)
  saveRDS(L$scale.i, here(vs.dir, "cov_scale.rds"))
  saveRDS(L$env.rct, here(vs.dir, "env_rct.rds"))
  saveRDS(L$env.rct.unscaled, here(vs.dir, "env_rct_unscaled.rds"))
  saveRDS(L$env.in, here(vs.dir, "env_in.rds"))
  saveRDS(L$env.args, here(vs.dir, "env_args.rds"))
  saveRDS(clim_X, here(vs.dir, "clim_X.rds"))
  saveRDS(p, here(vs.dir, "p.rds"))
  saveRDS(N_init, here(vs.dir, "N_init.rds"))
  saveRDS(sdd.pr, here(vs.dir, "sdd.rds"))
  saveRDS(sdd.ji, here(vs.dir, "sdd_ji.rds"))
  saveRDS(p.ji, here(vs.dir, "p_ji.rds"))
  saveRDS(U, here(vs.dir, "U.rds"))
  saveRDS(S, here(vs.dir, "S.rds"))
  saveRDS(lam.df, here(vs.dir, "lam_df.rds"))
}



