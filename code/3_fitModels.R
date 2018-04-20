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
sp <- "sp1"
overwrite <- TRUE
n_cores <- 4
issue <- c("none", "noise", "geogBias", "sampBias", 
           "noSB", "noDisp", "overDisp", "clim", "lc")[8]

# load workspace
pkgs <- c("dismo", "gbPopMod", "tidyverse", "magrittr", "MuMIn", "here", "doSNOW")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
p <- readRDS(here(paste0("out/", sp, "_p.rds")))
S <- readRDS(here(paste0("out/", sp, "_S.rds")))
U <- readRDS(here(paste0("out/", sp, "_U.rds")))
sdd.pr <- readRDS(here(paste0("out/", sp, "_sdd.rds")))
lam.df <- readRDS(here(paste0("out/", sp, "_lam_df.rds")))
env.in <- readRDS(here(paste0("out/", sp, "_env_in.rds")))
env.rct <- readRDS(here(paste0("out/", sp, "_env_rct.rds")))
env.rct.unsc <- readRDS(here(paste0("out/", sp, "_env_rct_unscaled.rds")))
n.cell <- nrow(env.in); n.grid <- nrow(env.rct)
issue_i <- read.csv(here("data/issues.csv"), stringsAsFactors=F)
sampling.issue <- filter(issue_i, Issue==issue)$Sampling
modeling.issue <- filter(issue_i, Issue==issue)$Modeling
O_Mx <- readRDS(here(paste0("out/", sp, "_O_Mx_", sampling.issue, ".rds")))
O_CA <- readRDS(here(paste0("out/", sp, "_O_CA_", sampling.issue, ".rds")))
O_IPM <- readRDS(here(paste0("out/", sp, "_O_IPM_", sampling.issue, ".rds")))


########
## Set model details
########
n_sim <- 8  # number of simulations per sample (mechanistic only)
v <- m <- n <- list(CA=NULL, IPM=NULL)

##--- CA
if(modeling.issue=="clim") {
  v$CA <- c("(Intercept)"=0, "temp"=0, "temp2"=0, "prec"=0, "prec2"=0)
} else if(modeling.issue=="lc") {
  v$CA <- c("(Intercept)"=0, "pOpn"=0, "pOth"=0, "pDec"=0, "pEvg"=0, "pMxd"=0)
} else {
  v$CA <- c("(Intercept)"=0, "temp"=0, "temp2"=0, "prec"=0, "prec2"=0, 
            "pOpn"=0, "pOth"=0, "pDec"=0)#, "pEvg"=0, "pMxd"=0)
}
m$CA <- paste(names(v$CA)[-1], collapse=" + ")
n$CA <- rep(list(length(v$CA)), 6)
names(n$CA) <- c("K", "s.jv", "s.ad", "p.f", "fec", "lam")
X.CA <- map(n$CA, ~cbind(1, as.matrix(env.in[,1:(.-1)])))

##--- IPM
if(modeling.issue=="clim") {
  v$IPM <- c("(Intercept)"=0, "size"=0, "size2"=0,
             "temp"=0, "temp2"=0, "prec"=0, "prec2"=0)
} else if(modeling.issue=="lc") {
  v$IPM <- c("(Intercept)"=0, "size"=0, "size2"=0, 
             "pOpn"=0, "pOth"=0, "pDec"=0, "pEvg"=0, "pMxd"=0)
} else {
  v$IPM <- c("(Intercept)"=0, "size"=0, "size2"=0,# "size3"=0, 
             "temp"=0, "temp2"=0, "prec"=0, "prec2"=0, 
             "pOpn"=0, "pOth"=0, "pDec"=0)#, "pEvg"=0, "pMxd"=0)
}
m$IPM <- paste(names(v$IPM)[-1], collapse=" + ")
n$IPM$z <- rep(list(sum(grepl("size", names(v$IPM))) + 1), 4)  # adds intercept
n$IPM$x <- rep(list(length(v$IPM) - n$IPM$z[[1]]), 4)
names(n$IPM$z) <- names(n$IPM$x) <- c("s", "g", "fl", "seed")
X.IPM <- map(n$IPM$x, ~as.matrix(env.in[,1:.]))


########
## fit models
########
##--- MaxEnt
cat("||||---- Beginning MaxEnt ---------------------------------------------\n")
Mx.f <- Mx.p <- vector("list", length(O_Mx))
for(i in 1:length(O_Mx)) {
  if(modeling.issue=="clim") {
    Mx.cov <- env.in[,c(1,3)]
  } else if(modeling.issue=="lc") {
    Mx.cov <- env.in[,5:9]
  } else {
    Mx.cov <- env.in[,c(1,3,5:9)]
  }
  Mx.f[[i]] <- maxent(x=Mx.cov, p=O_Mx[[i]])
  Mx.p[[i]] <- predict(Mx.f[[i]], Mx.cov)
}
P_Mx <- lam.df %>% select("x", "y", "x_y", "lat", "lon", "id", "id.inbd") %>%
  mutate(prP=apply(simplify2array(Mx.p), 1, mean)) %>%
  mutate(prP.sd=apply(simplify2array(Mx.p), 1, sd),
         Surv.S.f=prP/sum(prP)*sum(lam.df$Surv.S))

##--- CA
cat("||||---- Beginning CA -------------------------------------------------\n")
CA.f <- vector("list", length(O_CA))
for(i in 1:length(O_CA)) {
  O_CA.i <- O_CA[[i]]$d %>% filter(yr==p$tmax)
  O_CA.lam <- O_CA[[i]]$d %>% filter(!is.na(lambda))
  O_CA.K <- O_CA.i %>% filter(id %in% O_CA.lam$id[abs(O_CA.lam$lambda-1)<0.05])
  sim.ls <- sim.lam <- vector("list", n_sim)
  
  # global models
  options(na.action="na.fail")
  full.m <- list(K=glm(as.formula(paste("N ~", m$CA, collapse="")),
                       data=O_CA.K, family="poisson"),
                 s.jv=glm(as.formula(paste("cbind(s.jv.1, s.jv.0) ~", m$CA,
                                           collapse="")), 
                          data=O_CA.i, family="binomial"),
                 s.ad=glm(as.formula(paste("cbind(s.ad.1, s.ad.0) ~", m$CA,
                                           collapse="")), 
                          data=O_CA.i, family="binomial"),
                 p.f=glm(as.formula(paste("cbind(f.1, f.0) ~", m$CA,
                                          collapse="")), 
                         data=O_CA.i, family="binomial"),
                 fec=glm(as.formula(paste("fec ~", m$CA, collapse="")), 
                         data=O_CA.i, family="poisson"),
                 lam=lm(as.formula(paste("log(lambda) ~", m$CA, collapse="")),
                         data=O_CA.lam))
  
  # store coefficients from optimal models
  opt.m <- map(full.m, ~get.models(dredge(., m.lim=c(1,NA)), subset=1)[[1]])
  vars.opt <- map(opt.m, coef)
  vars.ls <- rep(list(v$CA), 6); names(vars.ls) <- names(opt.m)
  for(j in seq_along(vars.ls)) {
    vars.ls[[j]][names(vars.opt[[j]])] <- vars.opt[[j]]
  }
  
  # update parameters
  p.CA <- set_g_p(tmax=50, 
                  lc.r=diff(range(env.in$y)), lc.c=diff(range(env.in$x)),
                  n.lc=5, N.p.t0=n.cell, 
                  K=vars.ls$K, 
                  s.jv=vars.ls$s.jv,
                  s.ad=vars.ls$s.ad,
                  p.f=vars.ls$p.f,
                  fec=vars.ls$fec,
                  age.f=3, s.sb=p$s_SB, nSdFrt=1, 
                  p.est=as.matrix(logit(p$p_est)), 
                  sdd.max=p$sdd_max, sdd.rate=p$sdd_rate, n.ldd=1,
                  p.eat=as.matrix(1), bird.hab=p$bird_hab, s.bird=1, method="lm")
  p.CA$p_emig <- p$p_emig
  if(modeling.issue=="noSB") p.CA$s.sb <- 0
  if(modeling.issue=="noDisp") {
    p.CA$p_emig <- 0
    p.CA$sdd.rate <- 100
    p.CA$n.ldd <- 0
  }
  if(modeling.issue=="overDisp") {
    p.CA$p_emig <- min(5*p$p_emig, 0.99) 
    p.CA$sdd.rate <- p$sdd_rate/10
    p.CA$sdd.max <- p$sdd_max + 3
    p.CA$n.ldd <- 5
    sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, lc.col=8:12,
                            g.p=list(sdd.max=p.CA$sdd.max, 
                                     sdd.rate=p.CA$sdd.rate, 
                                     bird.hab=p.CA$bird.hab))
  }
  
  # run simulations
  CA.lc <- env.rct %>% rename(id.in=id.inbd)
  N.init <- matrix(0, n.grid, p.CA$age.f)  # column for each age class
  N.init[CA.lc$id[CA.lc$inbd], p.CA$age.f] <- p$n0
  N.init[CA.lc$id[CA.lc$inbd], -p.CA$age.f] <- round(p$n0/5)
  cat("||-- Starting simulations\n")
  if(n_cores > 1) {
    p.c <- makeCluster(n_cores); registerDoSNOW(p.c)
    sim.ls <- foreach(s=1:n_sim) %dopar% {
      gbPopMod::run_sim(n.grid, n.cell, p.CA, CA.lc, sdd.pr, N.init, NULL, F)
    }
    sim.lam <- foreach(s=1:n_sim) %dopar% {
      gbPopMod::run_sim_lambda(n.grid, n.cell, p.CA, vars.ls$lam, CA.lc, 
                     sdd.pr, N.init, "lm", F)
    }
    stopCluster(p.c)
  } else {
    for(s in 1:n_sim) {
      sim.ls[[s]] <- run_sim(n.grid, n.cell, p.CA, CA.lc, sdd.pr, N.init, NULL, F)
      sim.lam[[s]] <- run_sim_lambda(n.grid, n.cell, p.CA, vars.ls$lam, CA.lc, 
                                     sdd.pr, N.init, "lm", F)
    }
  }
  CA.f[[i]] <- aggregate_CA_simulations(sim.ls, p.CA$tmax, 
                                        max(p.CA$age.f), sim.lam)
  rm(sim.ls); rm(sim.lam)
  cat("  Finished dataset", i, "\n\n")
}
out <- summarize_CA_samples(CA.f, lam.df$id)
P_CAd <- lam.df %>% select("x", "y", "x_y", "lat", "lon", "id", "id.inbd") %>% 
  mutate(prP=out$prP[,p.CA$tmax+1],
         prP.sd=out$prP.sd[,p.CA$tmax+1],
         lam.S.f=rowMeans(out$N_ad.mn[,(-3:0)+p.CA$tmax]/
                            (out$N_ad.mn[,(-4:-1)+p.CA$tmax])),
         nSeed.f=out$nSd.mn[,p.CA$tmax], 
         D.f=out$D.mn[,p.CA$tmax],
         B0.f=out$B.mn[,1], 
         Btmax.f=out$B.mn[,p.CA$tmax+1],
         N.S.f=out$N_tot.mn[,p.CA$tmax+1], 
         Surv.S.f=out$N_ad.mn[,p.CA$tmax+1], 
         Rcr.S.f=out$N_rcr.mn[,p.CA$tmax+1],
         nSdStay.f=nSeed.f*(1-p.CA$p_emig), 
         nSdLeave.f=nSeed.f*p.CA$p_emig)
P_CAd$lam.S.f[is.nan(P_CAd$lam.S.f)] <- NA
P_CAl <- lam.df %>% select("x", "y", "x_y", "lat", "lon", "id", "id.inbd") %>% 
  mutate(prP=out$CA_lam.prP[,p.CA$tmax+1],
         prP.sd=out$CA_lam.prP.sd[,p.CA$tmax+1],
         Surv.S.f=out$CA_lam.N[,p.CA$tmax+1],
         lam.S.f=out$CA_lam.lam)

if(sum(is.na(P_CAd$Surv.S.f)>0)) cat("\n\n--------!! CA error\n\n")

##--- IPM
cat("||||---- Beginning IPM ------------------------------------------------\n")
set.seed(1)
S.f <- U.f <- vector("list", length(O_IPM))
for(i in 1:length(O_IPM)) {
  # separate data for MuMIn::dredge()
  O_IPM.i.s <- filter(O_IPM[[i]], !is.na(surv))
  O_IPM.i.g <- filter(O_IPM[[i]], !is.na(sizeNext) & !is.na(size))
  O_IPM.i.fl <- filter(O_IPM[[i]], !is.na(fl))
  O_IPM.i.seed <- filter(O_IPM[[i]], !is.na(seed))
  
  # full models
  options(na.action="na.fail")
  full.m <- list(s=glm(as.formula(paste("surv ~", m$IPM, collapse="")), 
                       data=O_IPM.i.s, family="binomial"),
                 g=lm(as.formula(paste("sizeNext ~", m$IPM, collapse="")), 
                      data=O_IPM.i.g),
                 fl=glm(as.formula(paste("fl ~", m$IPM, collapse="")), 
                        data=O_IPM.i.fl, family="binomial"),
                 seed=glm(as.formula(paste("seed ~", m$IPM, collapse="")), 
                          data=O_IPM.i.seed, family="poisson"))
  
  # store coefficients from optimal models
  opt.m <- map(full.m, ~get.models(dredge(., m.lim=c(1,NA)), subset=1)[[1]])
  vars.opt <- map(opt.m, coef)
  vars.ls <- rep(list(v$IPM), 4); names(vars.ls) <- names(opt.m)
  for(j in seq_along(vars.ls)) {
    vars.ls[[j]][names(vars.opt[[j]])] <- vars.opt[[j]]
  }
  
  # update parameters
  p.IPM <- p
  p.IPM$s_z <- vars.ls$s[1:n$IPM$z$s]
  p.IPM$s_x <- vars.ls$s[(n$IPM$z$s+1):length(v$IPM)]
  p.IPM$g_z <- vars.ls$g[1:n$IPM$z$g]
  p.IPM$g_x <- vars.ls$g[(n$IPM$z$g+1):length(v$IPM)]
  p.IPM$g_sig <- summary(opt.m$g)$sigma
  p.IPM$fl_z <- vars.ls$fl[1:n$IPM$z$fl]
  p.IPM$fl_x <- vars.ls$fl[(n$IPM$z$fl+1):length(v$IPM)]
  p.IPM$seed_z <- vars.ls$seed[1:n$IPM$z$seed]
  p.IPM$seed_x <- vars.ls$seed[(n$IPM$z$seed+1):length(v$IPM)]
  p.IPM$rcr_z <- filter(O_IPM[[i]], is.na(size)) %>% 
    summarise(mn=mean(sizeNext), sd=sd(sizeNext)) %>% unlist
  if(modeling.issue=="noSB") p.IPM$s_SB <- 0
  if(modeling.issue=="noDisp") {
    p.IPM$p_emig <- 0
    p.IPM$sdd.rate <- 100
    p.IPM$n.ldd <- 0
  }
  if(modeling.issue=="overDisp") {
    p.IPM$p_emig <- min(5*p$p_emig, 0.99) 
    p.IPM$sdd.rate <- p$sdd_rate/10
    p.IPM$sdd.max <- p$sdd_max + 3
    sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, lc.col=8:12,
                            g.p=list(sdd.max=p.IPM$sdd.max, 
                                     sdd.rate=p.IPM$sdd.rate, 
                                     bird.hab=p$bird_hab))
  }
  
  # use estimated slopes to fill IPM matrix
  cat("||-- Calculating IPM matrices\n")
  U.f[[i]] <- fill_IPM_matrices(n.cell, buffer=0.75, discrete=1, p.IPM, 
                                n$IPM$z, n$IPM$x, X.IPM, sdd.pr, env.in$id)
  
  # use estimated slopes to generate simulated data
  sim.ls <- vector("list", n_sim)
  cat("||-- Starting simulations\n")
  if(n_cores > 1) {
    p.c <- makeCluster(n_cores); registerDoSNOW(p.c)
    sim.ls <- foreach(s=1:n_sim) %dopar% {
      simulate_data(n.cell, U$lo, U$hi, p.IPM, X.IPM, n$IPM$z, sdd.pr, U$sdd.j)
    }
    stopCluster(p.c)
  } else {
    for(s in 1:n_sim) {
      sim.ls[[s]] <- simulate_data(n.cell, U$lo, U$hi, p.IPM, X.IPM, n$IPM$z, 
                                   sdd.pr, U$sdd.j)
      cat("||-- Finished simulation", s, "\n")
    }
  }
  S.f[[i]] <- aggregate_IPM_simulations(sim.ls, p.IPM$tmax)
  rm(sim.ls)
  cat("  Finished dataset", i, "\n\n")
}
out <- summarize_IPM_samples(U.f, S.f)

P_IPM <- lam.df %>% select("x", "y", "x_y", "lat", "lon", "id", "id.inbd") %>% 
  mutate(prP=out$Sf$prP,
         prP.sd=out$Sf$prP.sd,
         prP.U=out$Uf$prP[,p.IPM$tmax+1],
         lambda.f=apply(out$Uf$IPM.mn, 3, function(x) Re(eigen(x)$values[1])),
         lam.S.f=rowMeans(out$Sf$N_sim.mn[,(-3:0)+p.IPM$tmax]/
                            (out$Sf$N_sim.mn[,(-4:-1)+p.IPM$tmax]+0.0001)),
         nSeed.f=out$Sf$nSd.mn[,p.IPM$tmax], 
         D.f=out$Sf$D.mn[,p.IPM$tmax],
         B0.f=out$Sf$B.mn[,1], 
         Btmax.f=out$Sf$B.mn[,p.IPM$tmax+1],
         N.S.f=out$Sf$N_tot.mn, 
         Surv.S.f=out$Sf$N_surv.mn, 
         Rcr.S.f=out$Sf$N_rcr.mn,
         nSdStay.f=nSeed.f*(1-p.IPM$p_emig), 
         nSdLeave.f=nSeed.f*p.IPM$p_emig,
         N.U.f=apply(out$Uf$Nt.mn[,,p.IPM$tmax],2,sum), 
         lam.U.f=out$Uf$lam.mn[,p.IPM$tmax-1])

if(sum(is.na(P_IPM$Surv.S.f)>0)) cat("\n\n--------!! IPM error\n\n")

if(overwrite) {
  cat("Saving output\n")
  saveRDS(P_Mx, here(paste0("out/", sp, "_P_Mx_", issue, ".rds")))
  saveRDS(P_CAd, here(paste0("out/", sp, "_P_CAd_", issue, ".rds")))
  saveRDS(P_CAl, here(paste0("out/", sp, "_P_CAl_", issue, ".rds")))
  saveRDS(P_IPM, here(paste0("out/", sp, "_P_IPM_", issue, ".rds")))
}




