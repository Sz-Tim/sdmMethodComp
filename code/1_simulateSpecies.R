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
x_min <- 200
x_max <- Inf
y_min <- 0
y_max <- 85

# load workspace
pkgs <- c("gbPopMod", "tidyverse", "magrittr", "here", "doSNOW", "foreach")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
env.f <- paste0("data/ENF_", res, ".csv")
sp_i <- read.csv(paste0("data/species_", res, ".csv")) %>% filter(Name==sp)
vs.dir <- paste0("vs/", sp_i$Num)
nlcd.sp <- read.csv(here("data/PNAS_2017/", sp_i$LC_f))
L <- build_landscape(env.f, nlcd.sp, x_min, x_max, y_min, y_max) 
n.cell <- sum(L$env.rct$inbd)



########
## Set species traits
########
p <- fit_PNAS_species(sp, env.f, nlcd.sp, clim_X, FALSE, max_z_pow, habitat,
                      x_min, x_max, y_min, y_max)
if(sp=="garlic_mustard") {
  p$s_x <- c(-2, -0.9, -0.3, -0.6)
  p$g_x <- c(-4, -1.1, -1, -0.4)
  p$germ_x <- c(0.7, -0.5, -0.9, -0.2, -0.1)
  p$fl_x <- c(-0.3, -0.05, 0.25, -0.1)
  p$seed_x <- c(-1.3, -1, -0.3, -0.3)
  p$tmax <- 150
  p$bird_hab <- rep(1, 5)
  p$NDD_rcr <- 1  # mean number of recruits for NDD
  p$max_age <- 2
} else {
  p$s_x <- c(-5, -2.75, 1, -2.5)
  p$g_x <- c(-1.5, -0.4, -0.2, -0.5)
  p$germ_x <- c(-1.25, -4, -2, -2, -0.75)
  p$tmax <- 200
  p$bird_hab <- c(.32, .36, .05, .09, .09)
  p$NDD_rcr <- 40 # mean number of recruits for NDD
  p$max_age <- 50
}
p$n <- 30
p$tnonEq <- floor(p$tmax/3)
p$n0 <- 10
p$prop_init <- 0.0001
p$NDD <- T
p$sdd_max <- sp_i$sdd_max
p$sdd_rate <- sp_i$sdd_rate
p$ldd <- sp_i$ldd
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
X$germ <- cbind(1, X$germ[,-n_x$germ])

if(!all(file.exists(here(vs.dir, "sdd.rds")), 
        file.exists(here(vs.dir, "sdd_ji.rds")), 
        file.exists(here(vs.dir, "p_ji.rds")))) {
  cat("SDD neighborhoods not found. Calculating...\n")
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
} else {
  cat("SDD neighborhoods found\n")
  sdd.pr <- readRDS(here(vs.dir, "sdd.rds"))
  sdd.ji <- readRDS(here(vs.dir, "sdd_ji.rds"))
  p.ji <- readRDS(here(vs.dir, "p_ji.rds"))
}



########
## Generate underlying IPM and simulated data
########
# Initial populations
N_init <- rep(0, n.cell)
# N_init[sample(filter(L$env.in, x>425 & y>100 & y<200)$id.in,
#               p$prop_init*n.cell, replace=F)] <- p$n0
N_init[sample(filter(L$env.in, x>200 & y>40 & y<80)$id.in,
              p$prop_init*n.cell, replace=F)] <- p$n0

# Use assigned slopes to fill IPM matrix
U <- fill_IPM_matrices(n.cell, buffer=0, discrete=1, p, n_z, n_x, 
                       X, sdd.ji, p.ji, verbose=T)
if(sp=="garlic_mustard") {
  library(doSNOW); library(foreach)
  p.c <- makeCluster(n.cores); registerDoSNOW(p.c)
  U$lambda <- foreach(i=1:n.cell, .combine="c") %dopar% {
    iter_lambda(p, U$Ps[,,i], U$Fs[,,i], tol=0.5)
  }
  stopCluster(p.c)
} else {
  U$lambda <- apply(U$IPMs, 3, function(x) Re(eigen(x)$values[1]))
}
cat("Finished calculating U\n")

# Ground Truth: generate simulated data
S <- simulate_data(n.cell, U$lo, U$hi, p, X, n_z, sdd.ji, p.ji, N_init, sp, 
                   save_yrs=c(p$tnonEq+(-3:0), p$tmax+(-3:0)), T)
cat("Finished calculating S\n")

# Aggregate results
z.seq <- seq(p$z.rng[1], p$z.rng[2], length.out=5)
lam.df <- L$env.in %>%
  mutate(lambda=U$lambda,
         nSeed=S$nSd[,dim(S$nSd)[2]], 
         D=S$D[,dim(S$D)[2]], 
         B=S$B[,dim(S$B)[2]], 
         Surv.S=map_dbl(S$d, ~sum(.$surv[.$yr==p$tmax], na.rm=T)),
         Surv.S.t0=map_dbl(S$d, ~sum(.$surv[.$yr==p$tmax-1], na.rm=T)),
         Surv.lam=(Surv.S)/(Surv.S.t0),
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
         p_est.tmax=S$p_est.i[,dim(S$p_est.i)[2]],
         s.1=antilogit(c(cbind(1, z.seq[1], X$s) %*% c(p$s_z, p$s_x))),
         s.2=antilogit(c(cbind(1, z.seq[2], X$s) %*% c(p$s_z, p$s_x))),
         s.3=antilogit(c(cbind(1, z.seq[3], X$s) %*% c(p$s_z, p$s_x))),
         s.4=antilogit(c(cbind(1, z.seq[4], X$s) %*% c(p$s_z, p$s_x))),
         s.5=antilogit(c(cbind(1, z.seq[5], X$s) %*% c(p$s_z, p$s_x))),
         g=c(cbind(1, mean(p$z.rng), X$g) %*% c(p$g_z, p$g_x)),
         germ=antilogit(c(X$germ %*% p$germ_x)))
p$K_max <- max(lam.df$Surv.S)  # max abundance for CAd

library(viridis)
lam.gg <- ggplot(lam.df, aes(x=lon, y=lat)) + theme_bw() +
  theme(axis.text=element_blank()) + labs(x="", y="") +
  scale_fill_viridis(name="", option="B") +
  ggtitle(paste0(sp, ": 5km x 5km, favorable habitat"))

if(plots) {
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, lambda>1), aes(fill=lambda)) + 
    labs(subtitle="lambda") +
    geom_point(data=lam.df[which(N_init>0),], colour="white", shape=1)
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, Surv.S>0), aes(fill=Surv.S)) + 
    labs(subtitle=paste("N: year", p$tmax)) +
    geom_point(data=lam.df[which(N_init>0),], colour="white", shape=1)
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, Surv.S_nonEq>0), aes(fill=Surv.S_nonEq)) + 
    labs(subtitle=paste("N: year", p$tnonEq)) +
    geom_point(data=lam.df[which(N_init>0),], colour="white", shape=1)
  
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, nSeed>0), aes(fill=nSeed)) +
    labs(subtitle="Seed production") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, nRepro>0), aes(fill=nSeed/nRepro)) +
    labs(subtitle="Per capita seed production") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, D>0), aes(fill=round(D))) + 
    labs(subtitle="Immigrant seeds") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, D>0 | nSdStay>0), aes(fill=nSdStay+round(D))) + 
    labs(subtitle="Propagule pressure") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, B>0), aes(fill=B)) + 
    labs(subtitle="Seed bank") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(fill="gray50") +
    geom_tile(data=filter(lam.df, Rcr.S>0), aes(fill=Rcr.S)) +
    labs(subtitle="Recruits") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  
  ggplot() + theme_bw() +
    geom_tile(data=lam.df, aes(x=lon, y=lat), fill="gray90") +
    geom_tile(data=filter(lam.df, lambda>1), aes(x=lon, y=lat), 
              fill="darkblue", alpha=0.5) +
    geom_tile(data=filter(lam.df, Surv.S>0), aes(x=lon, y=lat), 
              fill="red", alpha=0.5) +
    theme(axis.text=element_blank()) + labs(x="", y="") +
    ggtitle(paste0(sp, ": 5km x 5km, favorable habitat"))
  
  lam.gg + geom_tile(aes(fill=s.1)) + 
    labs(subtitle="s.1") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=s.2)) + 
    labs(subtitle="s.2") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=s.3)) + 
    labs(subtitle="s.3") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=s.4)) + 
    labs(subtitle="s.4") +
    geom_point(data=lam.df[N_init>0,], colour="white", shape=1)
  lam.gg + geom_tile(aes(fill=s.5)) + 
    labs(subtitle="s.5") +
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



