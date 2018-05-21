# Functions for sampling, imposing issues, & fitting models
# Comparison of SDM approaches
# Tim Szewczyk



################################################################################
## Sampling functions
########

##-- CA sampling
##   For each sampled cell, arrange observed data for CA analysis
sample_for_CA <- function(S, Mech.sample, O_n, O_yr, n_samp, prop.sampled) {
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
    # p(establish) = recruits / attempted germinants
    # p_est = N_rcr / (nSeed*p_emig*rcr_dir + B*rcr_SB + D*rcr_dir)
    g <- with(O_CA[[s]], d$nSeed*p$p_emig*p$rcr_dir + B*p$rcr_SB + D*p$rcr_dir)
    O_CA[[s]]$d$p.est <- c(t(O_CA[[s]]$d$N.rcr / g))
  }
  return(O_CA)
}



##-- IPM sampling
##   For each sampled cell, arrange observed data for IPM analysis
sample_for_IPM <- function(S, Mech.sample, O_n, O_yr, n_samp, prop.sampled) {
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
  return(O_IPM)
}



##-- ISSUE::Noise
##   Add sampling error to observed datasets in the form of random noise
add_noise_Mx <- function(O_Mx, err, O_n, n_samp, Corr.sample, n.cell, P.i) {
  # Substitute false presences
  n.err <- err*O_n$Corr
  for(s in 1:n_samp) {
    Corr.sample[[s]][1:n.err] <- sample((1:n.cell)[-P.i], n.err, F)
  }
  return(map(Corr.sample, ~(1:n.cell %in% .)))
}
add_noise_CA <- function(O_CA, err, n_samp) {
  # N.obs = rnorm(N.true, N.true*N) | (N.obs ≥ 0)
  # fec.obs = rnorm(fec.true, fec.true*fec) | (fec.obs ≥ 0)
  for(s in 1:n_samp) {
    n_obs <- nrow(O_CA[[s]]$d)
    O_CA[[s]]$d <- O_CA[[s]]$d %>% 
      mutate(N=pmax(round(rnorm(n_obs, N, N*err$N)), 0),
             fec=pmax(round(rnorm(n_obs, fec, fec*err$fec))), 0) %>%
      group_by(id) %>% mutate(lambda=N/lag(N,1))
  }
  return(O_CA)
}
add_noise_IPM <- function(O_IPM, z_l, z_h, err, n_samp) {
  # sizeNext.obs = rnorm(SizeNext.true, g) | (z_l ≤ sizeNext.obs ≤ z_h)
  # seed.obs = rnorm(seed.true, seed.true*seed) | (seed.obs ≥ 0)
  options(warn=-1)
  for(s in 1:n_samp) {
    n_obs <- nrow(O_IPM[[s]])
    O_IPM[[s]] <- O_IPM[[s]] %>%
      mutate(sizeNext=pmin(rnorm(n_obs, sizeNext, err$g), z_h),
             sizeNext=pmax(sizeNext, z_l),
             seed=pmax(round(rnorm(n_obs, seed, seed*err$seed)), 0))
  }
  options(warn=0)
  return(O_IPM)
}




##-- ISSUE::Over Dispersal
##   Adjust dispersal parameters to overdisperse
add_misDisperse <- function(p.mod, p, sdd_max_adj=3, sdd_rate_adj=.1, ldd=5) {
  p.mod$sdd.max <- p.mod$sdd_max <- max(p$sdd_max + sdd_max_adj, 1)
  p.mod$sdd.rate <- p.mod$sdd_rate <- p$sdd_rate * sdd_rate_adj
  p.mod$p_emig <- pexp(0.5, p$sdd_rate, lower.tail=F)
  p.mod$n.ldd <- ldd
  return(p.mod)
}



################################################################################
## Fitting functions
########

##-- Fit MaxEnt
##
fit_Mx <- function(sp, sampling.issue, lam.df, vars, v.i) {
  library(dismo); library(here); library(tidyverse)
  # load observations
  O_Mx <- readRDS(here(paste0("out/", sp, "_O_Mx_", sampling.issue, ".rds")))
  X.Mx <- lam.df[,names(lam.df) %in% names(vars)[v.i]]
  # fit MaxEnt models
  Mx.f <- Mx.p <- vector("list", length(O_Mx))
  for(i in seq_along(O_Mx)) {
    Mx.f[[i]] <- maxent(x=X.Mx, p=O_Mx[[i]])
    Mx.p[[i]] <- predict(Mx.f[[i]], X.Mx)
  }
  # munge output
  P_Mx <- lam.df %>% select("x", "y", "x_y", "lat", "lon", "id", "id.inbd") %>%
    mutate(prP=apply(simplify2array(Mx.p), 1, mean)) %>%
    mutate(prP.sd=apply(simplify2array(Mx.p), 1, sd),
           Surv.S.f=prP/sum(prP)*sum(lam.df$Surv.S))
  return(P_Mx)
}



##--- Fit CA
##
fit_CA <- function(sp, sampling.issue, modeling.issue, p, env.rct, env.rct.unsc, 
                   lam.df, vars, v.i, v, m, N_init, sdd.pr, n.cell, n.grid, 
                   n_sim, n_cores) {
  library(here); library(tidyverse); library(gbPopMod); library(MuMIn); library(doSNOW)
  
  # load observations
  O_CA <- readRDS(here(paste0("out/", sp, "_O_CA_", sampling.issue, ".rds")))
  X.CA <- env.rct %>% rename(id.in=id.inbd) %>%
    select(one_of("x", "y", "x_y", "inbd", "id", "id.in", names(vars)[v.i]))
  
  # Fit CA models
  CA.f <- vector("list", length(O_CA))
  for(i in 1:length(O_CA)) {
    O_CA.i <- O_CA[[i]]$d %>% filter(yr==p$tmax)
    O_CA.lam <- O_CA[[i]]$d %>% filter(!is.na(lambda))
    O_CA.K <- O_CA.i %>% filter(id %in% O_CA.lam$id[abs(O_CA.lam$lambda-1)<0.05])
    sim.ls <- sim.lam <- vector("list", n_sim)
    
    # global models
    options(na.action="na.fail")
    full.m <- list(K=glm(as.formula(paste("N ~", m, collapse="")),
                         data=O_CA.K, family="poisson"),
                   s.jv=glm(as.formula(paste("cbind(s.jv.1, s.jv.0) ~", m,
                                             collapse="")), 
                            data=O_CA.i, family="binomial"),
                   s.ad=glm(as.formula(paste("cbind(s.ad.1, s.ad.0) ~", m,
                                             collapse="")), 
                            data=O_CA.i, family="binomial"),
                   p.f=glm(as.formula(paste("cbind(f.1, f.0) ~", m,
                                            collapse="")), 
                           data=O_CA.i, family="binomial"),
                   fec=glm(as.formula(paste("fec ~", m, collapse="")), 
                           data=O_CA.i, family="poisson"),
                   lam=lm(as.formula(paste("log(lambda) ~", m, collapse="")),
                          data=O_CA.lam))
    
    # store coefficients from optimal models
    opt.m <- map(full.m, ~get.models(dredge(., m.lim=c(1,NA)), subset=1)[[1]])
    vars.opt <- map(opt.m, coef)
    vars.ls <- rep(list(v), 6); names(vars.ls) <- names(opt.m)
    for(j in seq_along(vars.ls)) {
      vars.ls[[j]][names(vars.opt[[j]])] <- vars.opt[[j]]
    }
    
    # update parameters
    p.CA <- set_g_p(tmax=50, 
                    lc.r=diff(range(lam.df$y)), lc.c=diff(range(lam.df$x)),
                    n.lc=5, N.p.t0=n.cell, 
                    K=vars.ls$K, 
                    s.jv=vars.ls$s.jv,
                    s.ad=vars.ls$s.ad,
                    p.f=vars.ls$p.f,
                    fec=vars.ls$fec,
                    age.f=3, s.sb=p$s_SB, nSdFrt=1, 
                    p.est=as.matrix(logit(p$p_est)), 
                    sdd.max=p$sdd_max, sdd.rate=p$sdd_rate, n.ldd=1,
                    p.eat=matrix(1), bird.hab=p$bird_hab, s.bird=1, method="lm")
    p.CA$p_emig <- p$p_emig
    
    # impose issues
    if(modeling.issue=="noSB") p.CA$s.sb <- 0
    if(modeling.issue=="underDisp") {
      p.CA <- add_misDisperse(p.CA, p, sdd_max_adj=-2, sdd_rate_adj=10, ldd=1)
      sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, lc.col=8:12,
                              g.p=list(sdd.max=p.CA$sdd.max, 
                                       sdd.rate=p.CA$sdd.rate, 
                                       bird.hab=p.CA$bird.hab))
    }
    if(modeling.issue=="overDisp") {
      p.CA <- add_misDisperse(p.CA, p, sdd_max_adj=2, sdd_rate_adj=.1, ldd=3)
      sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, lc.col=8:12,
                              g.p=list(sdd.max=p.CA$sdd.max, 
                                       sdd.rate=p.CA$sdd.rate, 
                                       bird.hab=p.CA$bird.hab))
    }
    
    # run simulations
    N.init <- matrix(0, n.grid, p.CA$age.f)  # column for each age class
    N.init[X.CA$id[X.CA$inbd], p.CA$age.f] <- N_init
    N.init[X.CA$id[X.CA$inbd], -p.CA$age.f] <- round(N_init/5)
    cat("||-- Starting simulations\n")
    if(n_cores > 1) {
      p.c <- makeCluster(n_cores); registerDoSNOW(p.c)
      sim.ls <- foreach(s=1:n_sim) %dopar% {
        gbPopMod::run_sim(n.grid, n.cell, p.CA, X.CA, sdd.pr, N.init, NULL, F)
      }
      sim.lam <- foreach(s=1:n_sim) %dopar% {
        gbPopMod::run_sim_lambda(n.grid, n.cell, p.CA, vars.ls$lam, X.CA, 
                                 sdd.pr, N.init, "lm", F)
      }
      stopCluster(p.c)
    } else {
      for(s in 1:n_sim) {
        sim.ls[[s]] <- run_sim(n.grid, n.cell, p.CA, X.CA, sdd.pr, N.init, NULL, T)
        sim.lam[[s]] <- run_sim_lambda(n.grid, n.cell, p.CA, vars.ls$lam, X.CA, 
                                       sdd.pr, N.init, "lm", T)
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
  
  return(list(P_CAd=P_CAd, P_CAl=P_CAl))
}



##-- fit IPM
##
fit_IPM <- function(sp, sampling.issue, modeling.issue, p, env.rct.unsc, 
                   lam.df, vars, v.i, v, m, n_x, n_z, N_init, sdd.pr, n.cell,
                   n_sim, n_cores) {
  library(here); library(tidyverse); library(magrittr); library(MuMIn); library(doSNOW)
  walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), ~source(here(.)))
  
  # load observations
  O_IPM <- readRDS(here(paste0("out/", sp, "_O_IPM_", sampling.issue, ".rds")))
  X.IPM <- map(n_x, ~as.matrix(lam.df[,names(lam.df) %in% names(vars)[v.i]]))
  
  # Fit IPM models
  S.f <- U.f <- vector("list", length(O_IPM))
  for(i in 1:length(O_IPM)) {
    # separate data for MuMIn::dredge()
    O_IPM.i.s <- filter(O_IPM[[i]], !is.na(surv))
    O_IPM.i.g <- filter(O_IPM[[i]], !is.na(sizeNext) & !is.na(size))
    O_IPM.i.fl <- filter(O_IPM[[i]], !is.na(fl))
    O_IPM.i.seed <- filter(O_IPM[[i]], !is.na(seed))
    
    # full models
    options(na.action="na.fail")
    full.m <- list(s=glm(as.formula(paste("surv ~", m, collapse="")), 
                         data=O_IPM.i.s, family="binomial"),
                   g=lm(as.formula(paste("sizeNext ~", m, collapse="")), 
                        data=O_IPM.i.g),
                   fl=glm(as.formula(paste("fl ~", m, collapse="")), 
                          data=O_IPM.i.fl, family="binomial"),
                   seed=glm(as.formula(paste("seed ~", m, collapse="")), 
                            data=O_IPM.i.seed, family="poisson"))
    
    # store coefficients from optimal models
    opt.m <- map(full.m, ~get.models(dredge(., m.lim=c(1,NA)), subset=1)[[1]])
    vars.opt <- map(opt.m, coef)
    vars.ls <- rep(list(v), 4); names(vars.ls) <- names(opt.m)
    for(j in seq_along(vars.ls)) {
      vars.ls[[j]][names(vars.opt[[j]])] <- vars.opt[[j]]
    }
    
    # update parameters
    p.IPM <- p
    p.IPM$s_z <- vars.ls$s[1:n_z$s]
    p.IPM$s_x <- vars.ls$s[(n_z$s+1):length(v)]
    p.IPM$g_z <- vars.ls$g[1:n_z$g]
    p.IPM$g_x <- vars.ls$g[(n_z$g+1):length(v)]
    p.IPM$g_sig <- summary(opt.m$g)$sigma
    p.IPM$fl_z <- vars.ls$fl[1:n_z$fl]
    p.IPM$fl_x <- vars.ls$fl[(n_z$fl+1):length(v)]
    p.IPM$seed_z <- vars.ls$seed[1:n_z$seed]
    p.IPM$seed_x <- vars.ls$seed[(n_z$seed+1):length(v)]
    p.IPM$rcr_z <- filter(O_IPM[[i]], is.na(size)) %>% 
      summarise(mn=mean(sizeNext), sd=sd(sizeNext)) %>% unlist
    
    # impose issues
    if(modeling.issue=="noSB") p.IPM$s_SB <- 0
    if(modeling.issue=="underDisp") {
      p.IPM <- add_misDisperse(p.IPM, p, sdd_max_adj=-2, sdd_rate_adj=10, ldd=1)
      sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, lc.col=8:12,
                              g.p=list(sdd.max=p.IPM$sdd.max, 
                                       sdd.rate=p.IPM$sdd.rate, 
                                       bird.hab=p$bird_hab))
    }
    if(modeling.issue=="overDisp") {
      p.IPM <- add_misDisperse(p.IPM, p, sdd_max_adj=2, sdd_rate_adj=.1, ldd=5)
      sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, lc.col=8:12,
                              g.p=list(sdd.max=p.IPM$sdd.max, 
                                       sdd.rate=p.IPM$sdd.rate, 
                                       bird.hab=p$bird_hab))
    }
    
    # use estimated slopes to fill IPM matrix
    cat("||-- Calculating IPM matrices\n")
    U.f[[i]] <- fill_IPM_matrices(n.cell, buffer=0.75, discrete=1, p.IPM, n_z,
                                  n_x, X.IPM, sdd.pr, lam.df$id, N_init)
    
    # use estimated slopes to generate simulated data
    sim.ls <- vector("list", n_sim)
    cat("||-- Starting simulations\n")
    if(n_cores > 1) {
      p.c <- makeCluster(n_cores); registerDoSNOW(p.c)
      sim.ls <- foreach(s=1:n_sim) %dopar% {
        library(purrr); library(here)
        walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), ~source(here(.)))
        simulate_data(n.cell, U.f[[i]]$lo, U.f[[i]]$hi, p.IPM, X.IPM, n_z, 
                      sdd.pr, U.f[[i]]$sdd.j, N_init)
      }
      stopCluster(p.c)
    } else {
      for(s in 1:n_sim) {
        sim.ls[[s]] <- simulate_data(n.cell, U.f[[i]]$lo, U.f[[i]]$hi, p.IPM, X.IPM, n_z, 
                                     sdd.pr, U.f[[i]]$sdd.j, N_init)
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
  
  return(P_IPM)
}






