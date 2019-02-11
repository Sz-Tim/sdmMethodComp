# Functions for sampling, imposing issues, & fitting models
# Comparison of SDM approaches
# Tim Szewczyk



################################################################################
## Sampling functions
########

##-- sample for demographic CA model
#' Sample from S and summarise to population-level statistics for the specified
#' cells and years within each observed dataset
#' @param sp Virtual species name
#' @param S True individual-level data produced in 1_simulateSpecies.R
#' @param lam.df Dataframe with covariates, generated in 1_simulateSpecies.R
#' @param Mech.sample List of cell indices for each set of samples
#' @param O_yr List with [["CA"]] containing the years to sample
#' @param max_indiv Maximum number of individuals (per cell) to sample. A
#'   maximum of max_indiv/2 juveniles and max_indiv/2 adults are sampled
#' @return List with an element for each dataset, where each dataset is a list
#'   containing a dataframe "d" with population-level statistics for each cell
#'   and year, matrix "B" (dim=[n.cell, length(O_yr$CA)] with the size of the
#'   seedbank, and matrix "D" (dim=dim(B)) with the number of immigrant seeds to
#'   each cell
sample_for_CA <- function(sp, S, lam.df, Mech.sample, O_yr, max_indiv) {
  m <- ifelse(sp=="garlic_mustard", 2, 4)
  O_CA <- vector("list", length(Mech.sample))
  for(s in seq_along(O_CA)) {
    CA.d <- CA.B <- CA.D <- vector("list", length(Mech.sample[[1]]))
    germ.plots <- rbinom(length(Mech.sample[[s]]), 100, 
                         lam.df$germ[Mech.sample[[s]]])
    for(j in seq_along(Mech.sample[[s]])) {
      i <- Mech.sample[[s]][j]
      CA.d[[j]] <- data.frame(S$d[[i]]) %>% 
        filter(yr %in% O_yr$CA) %>%
        mutate(juv=age<m & age>0,
               adult=age>=m) %>%
        rbind(sample_n(filter(., juv), min(max_indiv$CA/2, sum(.$juv))),
              sample_n(filter(., adult), min(max_indiv$CA/2, sum(.$adult))))
      if(!all(O_yr$CA %in% CA.d[[j]]$yr)) {
        missing.yr <- O_yr$CA[!O_yr$CA %in% CA.d[[j]]$yr]
        CA.d[[j]] <- CA.d[[j]] %>%
          add_row(yr=missing.yr, 
                  age=NA, sizeNext=NA, seed=NA, fl=NA, surv=NA, size=NA)
      }
      CA.d[[j]] <- CA.d[[j]] %>% group_by(yr) %>% 
        summarise(N=sum(!is.na(size)), 
                  s.N.0=sum(surv[adult]==0, na.rm=TRUE),
                  s.N.1=sum(surv[adult]==1, na.rm=TRUE),
                  s.M.0=sum(surv[juv]==0, na.rm=TRUE),
                  s.M.1=sum(surv[juv]==1, na.rm=TRUE),
                  f.0=sum(fl[adult]==0, na.rm=TRUE),
                  f.1=sum(fl[adult]==1, na.rm=TRUE),
                  mu=median(seed, na.rm=TRUE) %>% round,
                  nSeed=sum(seed, na.rm=TRUE),
                  N.rcr=sum(is.na(size) & !is.na(sizeNext))) %>%
        mutate(mu=ifelse(is.nan(mu), 0, mu)) %>%
        ungroup() %>%
        mutate(lambda=N/lag(N,1)) %>%
        add_column(id.in=i) %>%
        full_join(lam.df[i,-(55:75)], by="id.in")
      CA.B[[j]] <- S$B[i,tail(1:dim(S$B)[2], length(O_yr$CA))]
      CA.D[[j]] <- S$D[i,tail(1:dim(S$B)[2], length(O_yr$CA))]
    }
    O_CA[[s]] <- list(d=do.call(rbind, CA.d),
                      B=do.call(rbind, CA.B),
                      D=do.call(rbind, CA.D))
    O_CA[[s]]$p_est <- lam.df[Mech.sample[[s]], -(55:75)]
    O_CA[[s]]$p_est$p.est.1 <- germ.plots
    O_CA[[s]]$p_est$p.est.0 <- 100-germ.plots
  }
  return(O_CA)
}



##-- sample for IPM and individual-based CA model
#' Sample from S for the specified cells and years within each observed dataset
#' @param p List of true parameters
#' @param S True individual-level data produced in 1_simulateSpecies.R
#' @param lam.df Dataframe with covariates, generated in 1_simulateSpecies.R
#' @param Mech.sample List of cell indices for each set of samples
#' @param O_yr List with [["IPM"]] containing the years to sample
#' @param max_indiv Maximum number of individuals (per cell) to sample. A
#'   maximum of max_indiv/3 each of small, medium, and large are sampled
#' @return List with an element for each dataset, where each dataset is a
#'   dataframe with individual-level statistics for each cell and year
sample_for_IPM <- function(p, S, lam.df, Mech.sample, O_yr, max_indiv) {
  O_IPM <- vector("list", length(Mech.sample))
  for(s in seq_along(O_IPM)) {
    IPM.d <- vector("list", length(Mech.sample[[1]]))
    germ.plots <- rbinom(length(Mech.sample[[s]]), 100, lam.df$germ[Mech.sample[[s]]])
    for(j in seq_along(Mech.sample[[s]])) {
      i <- Mech.sample[[s]][j]
      IPM.d[[j]] <- data.frame(S$d[[i]]) %>% filter(yr %in% O_yr$IPM) %>%
        mutate(size2=size^2, size3=size^3) %>%
        add_column(id.in=i) %>%
        full_join(env.in[i,], by="id.in")
      z_cut <- as.numeric(cut(IPM.d[[j]]$size, 
                              breaks=p$z.rng[1]+(0:3)*diff(p$z.rng)/3))
      sample.i <- c(which(is.na(IPM.d[[j]]$size)),
                    unlist(map(1:3, ~sample(which(z_cut==.), 
                                            min(max_indiv$IPM/3, 
                                                sum(z_cut==., na.rm=T))))))
      IPM.d[[j]] <- IPM.d[[j]][sample.i,]
    }
    O_IPM[[s]]$d <- do.call(rbind, IPM.d)
    O_IPM[[s]]$germ.df <- lam.df[Mech.sample[[s]],-(55:75)] %>%
      mutate(p.1=germ.plots, p.0=100-germ.plots)
  }
  return(O_IPM)
}






################################################################################
## Issue imposition functions
########

##-- NOISE: MaxEnt
#' Add sampling error to observed presence-absence datasets for MaxEnt, taking
#' the form of false presences
#' @param err Proportion of observations to make erroneous
#' @param Corr.sample List where each element is a vector of cell indexes
#' @param n.cell Number of inbound cells
#' @param P.i Cell indexes for true presences
#' @return List where each element is a logical vector (length n.cell)
#'   indicating observed presences and assumed absences after adding
#'   observation error
add_noise_Mx <- function(err, Corr.sample, n.cell, P.i) {
  n.err <- err*length(Corr.sample[[1]])
  for(s in seq_along(Corr.sample)) {
    Corr.sample[[s]][1:n.err] <- sample((1:n.cell)[-P.i], n.err, F)
  }
  return(map(Corr.sample, ~(1:n.cell %in% .)))
}



##-- NOISE: Demographic CA
#' Add sampling error to population-level metrics 
#' @param O_CA Output from \link{sample_for_CA()}
#' @param err List with elements [["N"]], describing the amount of error to add
#'   to observed abundance as N.obs = rnorm(N.true, N.true*N) subject to (N.obs
#'   ≥ 0), and [["mu"]], describing the amount of error to add to observed seed
#'   counts as mu.obs = rnorm(mu.true, mu.true*mu) subject to (mu.obs ≥ 0)
#' @return List of the same form as the output from \link{sample_for_CA()}
add_noise_CA <- function(O_CA, err) {
  for(s in seq_along(O_CA)) {
    n_obs <- nrow(O_CA[[s]]$d)
    O_CA[[s]]$d <- O_CA[[s]]$d %>% 
      mutate(N=pmax(round(rnorm(n_obs, N, N*err$N)), 0),
             mu=pmax(round(rnorm(n_obs, mu, mu*err$mu))), 0) %>%
      group_by(id) %>% mutate(lambda=N/lag(N,1))
  }
  return(O_CA)
}



##-- NOISE: Individual-based CA and IPM
#' Add sampling error to individual-level data
#' @param O_IPM Output from \link{sample_for_IPM()}
#' @param z_l Minimum allowable size
#' @param z_h Maximum allowable size
#' @param err List with elements [["g"]], describing the amount of error to add
#'   to observed growth as sizeNext.obs = rnorm(SizeNext.true, g) subject to
#'   (z_l ≤ sizeNext.obs ≤ z_h), and [["seed"]], describing the amount of error
#'   to add to observed seed counts as seed.obs = rnorm(seed.true,
#'   seed.true*seed) subject to (seed.obs ≥ 0)
#' @return List of the same form as the output from \link{sample_for_IPM()}
add_noise_IPM <- function(O_IPM, z_l, z_h, err) {
  options(warn=-1)
  for(s in seq_along(O_IPM)) {
    n_obs <- nrow(O_IPM[[s]]$d)
    O_IPM[[s]]$d <- O_IPM[[s]]$d %>%
      mutate(sizeNext=pmin(rnorm(n_obs, sizeNext, err$g), z_h),
             sizeNext=pmax(sizeNext, z_l),
             seed=pmax(round(rnorm(n_obs, seed, seed*err$seed)), 0))
  }
  options(warn=0)
  return(O_IPM)
}



##-- MISSPECIFIED DISPERSAL: Mechanistic models
#' Adjust dispersal parameters to over- or under-disperse
#' @param p.mod List of parameters fitted with data
#' @param p List of true parameters
#' @param sdd_max_adj Number of cells to adjust maximum distance for short
#'   distance dispersal; positive increases neighborhood size, negative reduces
#'   neighborhood size. Corrects to "1" when sdd_max is adjusted to < 1.
#' @param sdd_rate_adj Proportional adjustment to rate for short distance
#'   dispersal exponential kernel
#' @param ldd Number of annual long distance dispersal events
#' @return List with modified model parameters \code{p.mod}
add_misDisperse <- function(p.mod, p, sdd_max_adj=2, sdd_rate_adj=.1, ldd=5) {
  p.mod$sdd.max <- p.mod$sdd_max <- max(p$sdd_max + sdd_max_adj, 1)
  p.mod$sdd.rate <- p.mod$sdd_rate <- p$sdd_rate * sdd_rate_adj
  p.mod$p_emig <- pexp(0.5, p.mod$sdd_rate, lower.tail=F)
  p.mod$n.ldd <- p.mod$ldd <- ldd
  return(p.mod)
}






################################################################################
## Fitting functions
########

##-- fit MaxEnt
#' Generate predicted distribution using MaxEnt
#' @param spNum Virtual species number
#' @param issue Issue name
#' @param samp.issue Issue with observed dataset
#' @param lam.df Dataframe with covariates, generated in 1_simulateSpecies.R
#' @param v Variables to include in model
#' @return List with dataframe [["P_MxE"]] containing a row for each cell and
#'   columns for cell information, mean(probability of presence), and
#'   sd(probability of presence); and [["diag"]] containing diagnostics
fit_MxE <- function(spNum, issue, samp.issue, lam.df, v) {
  library(here); library(tidyverse)
  path_iss <- paste0("out/maxent/", spNum, "/", issue, "/")
  # load observations
  O_Mx <- readRDS(here("vs", spNum, paste0("O_Mx_", samp.issue, ".rds")))
  X.Mx <- lam.df[,names(lam.df) %in% v]
  temp.rast <- vector("list", ncol(X.Mx))
  for(vi in 1:ncol(X.Mx)) {
    temp.rast[[vi]] <- raster::rasterFromXYZ(cbind(lam.df[,c("lon","lat")], X.Mx[,vi]))
  }
  names(temp.rast) <- names(X.Mx)
  rast.Mx <- raster::stack(temp.rast)

  bias_type="none"
  MxE.f <- MxE.p <- vector("list", length(O_Mx))
  fit.args <- c("responsecurves", "jackknife", "replicates=10", "plots",
                "removeduplicates", "nothreshold", "nohinge", "pictures",
                "noautofeature", "noprefixes", "writeplotdata",
                "outputformat=logistic")
  for(i in seq_along(O_Mx)) {
    if(!dir.exists(paste0(path_iss, i))) dir.create(paste0(path_iss, i), recursive=T)
    MxE.f[[i]] <- dismo::maxent(x=rast.Mx, 
                                p=as.matrix(lam.df[O_Mx[[i]], c("lon", "lat")]),
                                args=fit.args, path=paste0(path_iss, i))
    MxE.p[[i]] <- raster::calc(dismo::predict(MxE.f[[i]], rast.Mx, 
                                              args="outputformat=logistic"), 
                               fun=mean, na.rm=T)
    thresh_j <- rep(NA, 10)
    for (j in 0:9){   #loop for each replicate
      d <- read.csv(paste0(path_iss, i, "/species_", j, "_samplePredictions.csv"))
      thresh_j[j+1] <- min(dplyr::filter(d, Test.or.train=="train")$Logistic.prediction)
    }
    MxE.p[[i]] <- raster::setValues(MxE.p[[i]],
                                    raster::values(MxE.p[[i]]) > mean(thresh_j))
  }
  # munge output
  names(MxE.p) <- 1:length(MxE.p)
  MxE.p.data <- map_dfr(MxE.p, ~na.omit(.@data@values)) %>% as.matrix
  P_MxE <- lam.df %>% 
    dplyr::select("x", "y", "x_y", "lat", "lon", "id", "id.in") %>%
    mutate(prP=apply(MxE.p.data, 1, mean),
           prP.sd=apply(MxE.p.data, 1, sd))
  diagnostics <- NULL
  TSS_MxE <- list(N=apply(MxE.p.data, 2, calc_TSS, S.pa=lam.df$Surv.S>0),
                  lam=apply(MxE.p.data, 2, calc_TSS, S.pa=lam.df$lambda>1))
  # diagnostics <- map(MxE.f, ~evaluate(p=S_p, a=S_a, model=., x=rast.Mx))
  return(list(P_MxE=P_MxE, diag=diagnostics, TSS_MxE=TSS_MxE))
}



##-- fit MaxLike
#' Generate predicted distribution using MaxLike
#' @param sp Virtual species
#' @param issue Issue name
#' @param samp.issue Issue with observed dataset
#' @param lam.df Dataframe with covariates, generated in 1_simulateSpecies.R
#' @param v Variables to include in model
#' @param m Model formula
#' @return List with dataframe [["P_MxL"]] containing a row for each cell and
#'   columns for cell information, thresholded mean(probability of presence),
#'   raw mean(probability of presence), and sd(probability of presence); and
#'   [["diag"]] containing diagnostics
fit_MxL <- function(sp, issue, samp.issue, lam.df, v, m) {
  library(here); library(tidyverse); library(raster); library(maxlike)
  
  O_Mx <- readRDS(here("vs/", sp, paste0("O_Mx_", samp.issue, ".rds")))
  X.Mx <- lam.df[,names(lam.df) %in% v]
  temp.rast <- vector("list", ncol(X.Mx))
  for(vi in 1:ncol(X.Mx)) {
    temp.rast[[vi]] <- rasterFromXYZ(cbind(lam.df[,c("lon", "lat")], X.Mx[,vi]))
  }
  names(temp.rast) <- names(X.Mx)
  rast.Mx <- stack(temp.rast)
  MxL.f <- MxL.p <- MxL.PA <- thresh <- vector("list", length(O_Mx))
  for(i in  seq_along(O_Mx)) {
    MxL.f[[i]] <- maxlike(as.formula(paste("~", m)), rast.Mx, hessian=F,
                         as.matrix(lam.df[O_Mx[[i]], c("lon", "lat")]), savedata=T)
    MxL.p[[i]] <- predict(MxL.f[[i]])
    thresh[[i]] <- min(MxL.p[[i]]@data@values[lam.df$id[O_Mx[[i]]]])
    MxL.PA[[i]] <- MxL.p[[i]]@data@values > thresh[[i]]
  }
  names(MxL.p) <- 1:length(MxL.p)
  MxL.p.data <- map_dfr(MxL.p, ~.@data@values) %>% as.matrix
  MxL.PA <- do.call("cbind", MxL.PA)
  P_MxL <- lam.df %>% 
    dplyr::select("x", "y", "x_y", "lat", "lon", "id", "id.in") %>%
    mutate(prP=apply(MxL.PA[lam.df$id,], 1, mean),
           prP_raw=apply(MxL.p.data[lam.df$id,], 1, mean),
           prP.sd=apply(MxL.p.data[lam.df$id,], 1, sd))
  diagnostics <- NULL#map(MxL.f, summary)
  return(list(P_MxL=P_MxL, diag=diagnostics))
}



##-- fit demographic CA model
#' Generate predicted distribution using a demographic CA model
#' @param sp Virtual species name
#' @param sp_i Single row dataframe with species information
#' @param samp.issue Issue with observed dataset
#' @param mod.issue Issue to impose on model
#' @param p List of true parameters
#' @param env.rct Dataframe with scaled covariates for full rectangular grid
#' @param env.rct.unsc Dataframe with raw covariates for full rectangular grid
#' @param lam.df Dataframe with covariates, generated in 1_simulateSpecies.R
#' @param v Variables to include in full model
#' @param m Model formula for full model in each population-level regression
#' @param N_init Vector with initial population size in each cell
#' @param sdd.pr List with short distance dispersal neighborhoods generated in
#'   1_simulateSpecies.R
#' @param sdd.ji
#' @param p.ji
#' @param n_sim Number of simulations to run per observed dataset
#' @param n_cores Number of cores to use for running simulations in parallel
#' @return List with dataframe P_CAi with overall summaries across datasets for
#'   the demographic CA distribution and diagnostics containing the vital rate
#'   regressions
fit_CA <- function(sp, sp_i, samp.issue, mod.issue, p, env.rct, env.rct.unsc, 
                   lam.df, v, m, N_init, sdd.pr, sdd.ji, p.ji, n_sim, n_cores) {
  library(here); library(tidyverse); library(gbPopMod); 
  library(MuMIn); library(lme4); library(doSNOW)
  walk(dir("code", "fn", full.names=T), source)
  out.dir <- paste0("out/", sp_i$Num, "/", samp.issue, "/", mod.issue, "/")
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  n.cell <- nrow(lam.df)
  n.grid <- nrow(env.rct)
  
  # load observations
  O_CA <- readRDS(here("vs", sp_i$Num, paste0("O_CA_", samp.issue, ".rds")))
  X.CA <- env.rct %>% rename(id.in=id.in) %>%
    dplyr::select(one_of("x", "y", "x_y", "inbd", "id", "id.in", names(v)[-1]))
  n_LC <- n_distinct(read.csv(paste0("data/PNAS_2017/", sp_i$LC_f))$agg)
  
  # initialize CA parameters
  p.CA <- set_g_p(tmax=p$tmax, 
                  lc.r=diff(range(lam.df$y)), lc.c=diff(range(lam.df$x)),
                  n.lc=5, N.p.t0=n.cell, sdd.st=F,
                  m=sp_i$m, gamma=1, 
                  s.B=p$s_SB, g.D=p$rcr_dir, g.B=p$rcr_SB,
                  sdd.max=p$sdd_max, sdd.rate=p$sdd_rate, n.ldd=p$ldd,
                  p.c=matrix(1), bird.hab=p$bird_hab, s.c=0.6, method="lm")
  p.CA$p_emig <- p$p_emig
  
  # impose dispersal issues
  if(mod.issue=="underDisp") {
    p.CA <- add_misDisperse(p.CA, p, sdd_max_adj=-2, sdd_rate_adj=2, 
                            ldd=round(p$ldd/5))
  } else if(mod.issue=="overDisp") {
    p.CA <- add_misDisperse(p.CA, p, sdd_max_adj=2, sdd_rate_adj=.5, 
                            ldd=p$ldd*5)
  }
  if(grepl("Disp", mod.issue)) {
    sdd.pr <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, 
                            lc.col=tail(1:ncol(env.rct.unsc), n_LC),
                            g.p=list(sdd.max=p.CA$sdd.max, 
                                     sdd.rate=p.CA$sdd.rate, 
                                     bird.hab=p.CA$bird.hab))
    sdd.ji.rows <- lapply(1:n.cell, function(x) which(sdd.pr$sp.df$j.idin==x))
    sdd.ji <- lapply(sdd.ji.rows, function(x) sdd.pr$sp.df$i.idin[x]) 
    p.ji <- lapply(sdd.ji.rows, function(x) sdd.pr$sp.df$pr[x]) 
  }
  
  # Fit CA models
  p.c <- makeCluster(n_cores); registerDoSNOW(p.c)
  foreach(i=1:length(O_CA), .errorhandling="pass", 
          .packages=c("tidyverse", "gbPopMod", "MuMIn", "lme4")) %dopar% {
    walk(dir("code", "fn", full.names=T), source)
    i_pad <- str_pad(i, 2, pad="0")
    O_CA.i <- O_CA[[i]]$d 
    O_CA.p_est <- O_CA[[i]]$p_est
    O_CA.K <- O_CA.i %>% filter(id %in% O_CA.i$id[abs(O_CA.i$lambda-1)<0.05])
    sim.ls <- vector("list", n_sim)
    
    # global models
    options(na.action="na.fail")
    full.m <- opt.m <- vars.opt <- setNames(vector("list", 6), 
                                            c("K", "s.M", "s.N", "p.f", "mu", "p.est"))
    
    full.m$K <- glmer(as.formula(paste("N ~", m, collapse="")),
                      data=O_CA.K, family="poisson")
    if(sum(O_CA.i$s.M.0==0)/nrow(O_CA.i) < 0.9) { # in case very low mortality
      full.m$s.M <- glmer(as.formula(paste("cbind(s.M.1, s.M.0) ~", m, collapse="")), 
                          data=O_CA.i, family="binomial")
    }
    if(sp=="barberry" && sum(O_CA.i$s.N.0==0)/nrow(O_CA.i) < 0.9) { # <10% mortality
      full.m$s.N <- glmer(as.formula(paste("cbind(s.N.1, s.N.0) ~", m, collapse="")), 
                          data=O_CA.i, family="binomial")
    } 
    full.m$p.f <- glmer(as.formula(paste("cbind(f.1, f.0) ~", m, collapse="")), 
                        data=O_CA.i, family="binomial")
    full.m$mu <- glmer(as.formula(paste("mu ~", m, collapse="")), 
                       data=O_CA.i, family="poisson")
    full.m$p.est <- glm(as.formula(paste("cbind(p.est.1, p.est.0) ~", 
                                         substr(m, 1, nchar(m)-9), 
                                         collapse="")), 
                        data=O_CA.p_est, family="binomial")
    
    # store coefficients from optimal models
    vars.ls <- rep(list(v), length(full.m)); names(vars.ls) <- names(opt.m)
    for(j in seq_along(full.m)) {
      if(is.null(full.m[[j]])) {
        if(sp=="barberry") {
          vars.opt[[j]] <- c("(Intercept)"=logit(0.9999))
        } else if(sp=="mustard") {
          vars.opt[[j]] <- c("(Intercept)"=logit(0.0001))
        }
      } else {
        opt.m[[j]] <- get.models(dredge(full.m[[j]], m.max=6), subset=1)[[1]]
        if(class(coef(opt.m[[j]]))=="numeric") {
          vars.opt[[j]] <- coef(opt.m[[j]])
        } else {
          vars.opt[[j]] <- colMeans(coef(opt.m[[j]])$yr)
        }
      }
      vars.ls[[j]][names(vars.opt[[j]])] <- vars.opt[[j]]
    }
    
    # update parameters
    p.CA$K <- vars.ls$K
    p.CA$s.M <- vars.ls$s.M
    p.CA$s.N <- vars.ls$s.N
    p.CA$p.f <- vars.ls$p.f
    p.CA$mu <- vars.ls$mu
    p.CA$p <- vars.ls$p.est
    
    # impose seed bank issue
    if(mod.issue=="noSB") {p.CA$s.B <- 0; p.CA$bank=F}
    
    # save parameters to diagnostic file
    saveRDS(p.CA, paste0(out.dir, "/CAd_diag_", i_pad, ".rds"))
    
    # run simulations
    N.init <- matrix(0, n.grid, p.CA$m)  # column for each age class
    N.init[X.CA$id[X.CA$inbd], p.CA$m] <- N_init
    N.init[X.CA$id[X.CA$inbd], -p.CA$m] <- round(N_init/5)
    for(s in 1:n_sim) {
      sim.ls[[s]] <- gbPopMod::run_sim(n.grid, n.cell, p.CA, X.CA, sdd.pr, 
                                       N.init, NULL, T, (-1:0)+p.CA$tmax, 
                                       p$K_max, dem_out=TRUE, FALSE)
    }
    saveRDS(aggregate_CAd_simulations(sim.ls, max(p.CA$m)), 
            paste0(out.dir, "/CAd_fit_", i_pad, ".rds"))
  }
  stopCluster(p.c)
  
  fits <- list.files(out.dir, "CAd_fit", full.names=T) %>% map(readRDS) 
  PA.mx <- do.call("cbind", map(fits, ~.$P[,dim(fits[[1]]$P)[2],1]))[lam.df$id,]
  out <- summarize_CAd_samples(fits, lam.df$id)
  TSS_CAd <- list(N=apply(PA.mx, 2, calc_TSS, S.pa=lam.df$Surv.S>0),
                  lam=apply(PA.mx, 2, calc_TSS, S.pa=lam.df$lambda>1))
  diagnostics <- list.files(out.dir, "CAd_diag", full.names=T) %>% map(readRDS)
  P_CAd <- lam.df %>% 
    dplyr::select("x", "y", "x_y", "lat", "lon", "id", "id.in") %>% 
    mutate(prP=out$prP[,dim(out$prP)[2]],
           nSeed.f=out$nSd.mn[,dim(out$nSd.mn)[2]], 
           D.f=out$D.mn[,dim(out$D.mn)[2]],
           B.f=out$B.mn[,dim(out$B.mn)[2]],
           N.S.f=out$N_tot.mn[,dim(out$N_tot.mn)[2]], 
           Surv.S.f=out$N_ad.mn[,dim(out$N_ad.mn)[2]], 
           Rcr.S.f=out$N_rcr.mn[,dim(out$N_rcr.mn)[2]],
           nSdStay.f=nSeed.f*(1-p.CA$p_emig), 
           nSdLeave.f=nSeed.f*p.CA$p_emig)
  
  return(list(P_CAd=P_CAd, diag=diagnostics, TSS_CAd=TSS_CAd))
}



##-- fit IPM & individual-based CA model
#' Generate predicted distribution using an IPM and an individual-based CA model
#' @param sp Virtual species name
#' @param sp_i Single row dataframe with species information
#' @param samp.issue Issue with observed dataset
#' @param mod.issue Issue to impose on model
#' @param p List of true parameters
#' @param env.rct.unsc Dataframe with raw covariates for full rectangular grid
#' @param lam.df Dataframe with covariates, generated in 1_simulateSpecies.R
#' @param v Variables to include in full model
#' @param m Model formula for full model in each population-level regression
#' @param n_z List of number of size covariates for each vital rate regression
#' @param n_x List of number of environmental covariates for each vital rate 
#' regression
#' @param N_init Vector with initial population size in each cell
#' @param sdd.ji
#' @param p.ji
#' @param n_sim Number of simulations to run per observed dataset
#' @param n_cores Number of cores to use for running simulations in parallel
#' @return List with dataframes P_IPM and P_CAi with overall summaries across
#'   datasets for the IPM distribution and individual-based CA distribution
#'   respectively, and diagnostics containing the vital rate regressions
fit_IPM <- function(sp, sp_i, samp.issue, mod.issue, p, env.rct.unsc, lam.df, v,
                    m, n_z, n_x, N_init, sdd.ji, p.ji, n_sim, n_cores) {
  library(here); library(tidyverse); library(magrittr); library(gbPopMod);
  library(MuMIn); library(doSNOW)
  walk(dir("code", "fn", full.names=T), source)
  out.dir <- paste0("out/", sp_i$Num, "/", samp.issue, "/", mod.issue, "/")
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  n.cell <- nrow(lam.df)
  
  # load observations
  O_IPM <- readRDS(here("vs", sp_i$Num, paste0("O_IPM_", samp.issue, ".rds")))
  X.IPM <- map(n_x, ~as.matrix(lam.df[,names(lam.df) %in% names(v)[-(1:n_z$s)]]))
  if(!is.null(X.IPM$germ)) X.IPM$germ <- cbind(1, X.IPM$germ)
  
  # initialize IPM parameters
  p.IPM <- p
  
  # impose dispersal issues
  if(mod.issue=="underDisp") {
    p.IPM <- add_misDisperse(p.IPM, p, sdd_max_adj=-2, sdd_rate_adj=2, 
                             ldd=round(p$ldd/5))
  } else if(mod.issue=="overDisp") {
    p.IPM <- add_misDisperse(p.IPM, p, sdd_max_adj=2, sdd_rate_adj=.5, 
                             ldd=p$ldd*5)
  }
  if(grepl("Disp", mod.issue)) {
    sp.df <- sdd_set_probs(ncell=n.cell, lc.df=env.rct.unsc, lc.col=8:12,
                           g.p=list(sdd.max=p.IPM$sdd.max, 
                                    sdd.rate=p.IPM$sdd.rate, 
                                    bird.hab=p$bird_hab))$sp.df
    sdd.ji.rows <- lapply(1:n.cell, function(x) which(sp.df$j.idin==x))
    sdd.ji <- lapply(sdd.ji.rows, function(x) sp.df$i.idin[x]) 
    p.ji <- lapply(sdd.ji.rows, function(x) sp.df$pr[x]) 
  }
  
  # Fit IPM models
  p.c <- makeCluster(n_cores); registerDoSNOW(p.c)
  foreach(i=1:length(O_IPM), .errorhandling="pass", 
          .packages=c("tidyverse", "magrittr", "gbPopMod", "MuMIn")) %dopar% {
    walk(dir("code", "fn", full.names=T), source)
    sim.ls <- vector("list", n_sim)
    i_pad <- str_pad(i, 2, pad="0")
    # separate data for MuMIn::dredge()
    O_IPM.i.s <- filter(O_IPM[[i]]$d, !is.na(surv))
    if(sp=="garlic_mustard") O_IPM.i.s <- filter(O_IPM.i.s, age==1)
    O_IPM.i.g <- filter(O_IPM[[i]]$d, !is.na(sizeNext) & !is.na(size))
    O_IPM.i.fl <- filter(O_IPM[[i]]$d, !is.na(fl))
    O_IPM.i.seed <- filter(O_IPM[[i]]$d, !is.na(seed))
    O_IPM.i.germ <- O_IPM[[i]]$germ.df
    
    # full models
    options(na.action="na.fail")
    full.m <- list(s=glm(as.formula(paste("surv ~", m, collapse="")), 
                         data=O_IPM.i.s, family="binomial"),
                   g=lm(as.formula(paste("sizeNext ~", m, collapse="")), 
                        data=O_IPM.i.g),
                   fl=glm(as.formula(paste("fl ~", m, collapse="")), 
                          data=O_IPM.i.fl, family="binomial"),
                   seed=glm(as.formula(paste("seed ~", m, collapse="")), 
                            data=O_IPM.i.seed, family="poisson"),
                   germ=glm(as.formula(paste("cbind(p.1, p.0) ~", 
                                             paste(grep("size", names(v), 
                                                        invert=T, value=T)[-1], 
                                                   collapse=" + "), 
                                             collapse="")),
                            data=O_IPM.i.germ, family="binomial"))
    
    # store coefficients from optimal models
    opt.m <- map(full.m, ~get.models(dredge(., m.lim=c(1,6)), subset=1)[[1]])
    vars.opt <- map(opt.m, coef)
    vars.ls <- rep(list(v), 5); names(vars.ls) <- names(opt.m)
    for(j in seq_along(vars.ls)) {
      vars.ls[[j]][names(vars.opt[[j]])] <- vars.opt[[j]]
    }
    
    # update parameters
    p.IPM$s_z <- vars.ls$s[1:n_z$s]
    p.IPM$s_x <- vars.ls$s[(n_z$s+1):length(v)]
    p.IPM$g_z <- vars.ls$g[1:n_z$g]
    p.IPM$g_x <- vars.ls$g[(n_z$g+1):length(v)]
    p.IPM$g_sig <- summary(opt.m$g)$sigma
    p.IPM$fl_z <- vars.ls$fl[1:n_z$fl]
    p.IPM$fl_x <- vars.ls$fl[(n_z$fl+1):length(v)]
    p.IPM$seed_z <- vars.ls$seed[1:n_z$seed]
    p.IPM$seed_x <- vars.ls$seed[(n_z$seed+1):length(v)]
    p.IPM$germ_x <- vars.ls$germ[grep("size", names(v), invert=T, value=T)]
    p.IPM$rcr_z <- filter(O_IPM[[i]]$d, is.na(size)) %>% 
      summarise(mn=mean(sizeNext), sd=sd(sizeNext)) %>% unlist
    
    # impose issues
    if(mod.issue=="noSB") p.IPM$s_SB <- 0
    
    # save parameters to diagnostics file
    saveRDS(p.IPM, paste0(out.dir, "/IPM_diag_", i_pad, ".rds"))
    
    # use estimated slopes to fill IPM matrix
    U.f <- fill_IPM_matrices(n.cell, buffer=0, discrete=1, p.IPM, n_z,
                                  n_x, X.IPM, sdd.ji, p.ji, sp)
    if(sp=="garlic_mustard") {
      U.f$lambda <- sapply(1:n.cell, 
                         function(x) iter_lambda(p.IPM, U.f$Ps[,,x], U.f$Fs[,,x]))
    } else {
      U.f$lambda <- apply(U.f$IPMs, 3, function(x) Re(eigen(x)$values[1]))
    }
    saveRDS(U.f, paste0(out.dir, "/IPM_fit_", i_pad, ".rds"))
    
    # use estimated slopes to generate simulated data
    for(s in 1:n_sim) {
      sim.ls[[s]] <- simulate_data(n.cell, U.f$lo, U.f$hi, p.IPM, X.IPM, n_z, 
                                   sdd.ji, p.ji, N_init, sp, save_yrs=p.IPM$tmax)
    }
    saveRDS(aggregate_CAi_simulations(sim.ls, p.IPM$tmax), 
            paste0(out.dir, "/CAi_fit_", i_pad, ".rds"))
  }
  stopCluster(p.c)
  
  out <- summarize_IPM_CAi_samples(
    map(list.files(out.dir, "IPM_fit", full.names=T), readRDS),
    map(list.files(out.dir, "CAi_fit", full.names=T), readRDS)
  )
  TSS_IPM <- list(N=apply(out$Uf.pa, 2, calc_TSS, S.pa=lam.df$Surv.S>0),
                  lam=apply(out$Uf.pa, 2, calc_TSS, S.pa=lam.df$lambda>1))
  TSS_CAi <- list(N=apply(out$Sf.pa, 2, calc_TSS, S.pa=lam.df$Surv.S>0),
                  lam=apply(out$Sf.pa, 2, calc_TSS, S.pa=lam.df$lambda>1))
  diagnostics <- list.files(out.dir, "IPM_diag", full.names=T) %>% map(readRDS)
  
  P_CAi <- lam.df %>% 
    dplyr::select("x", "y", "x_y", "lat", "lon", "id", "id.in") %>% 
    mutate(prP=out$Sf$prP,
           nSeed.f=out$Sf$nSd.mn[,dim(out$Sf$nSd.mn)[2]], 
           D.f=out$Sf$D.mn[,dim(out$Sf$D.mn)[2]],
           B.f=out$Sf$B.mn[,dim(out$Sf$B.mn)[2]],
           N.S.f=out$Sf$N_tot.mn, 
           Surv.S.f=out$Sf$N_surv.mn, 
           Rcr.S.f=out$Sf$N_rcr.mn,
           nSdStay.f=nSeed.f*(1-p.IPM$p_emig), 
           nSdLeave.f=nSeed.f*p.IPM$p_emig)
  P_IPM <- lam.df %>% 
    dplyr::select("x", "y", "x_y", "lat", "lon", "id", "id.in") %>% 
    mutate(prP=out$Uf$prP,
           lambda.f=out$Uf$lam.mn)
  
  return(list(P_IPM=P_IPM, P_CAi=P_CAi, diag=diagnostics, 
              TSS_IPM=TSS_IPM, TSS_CAi=TSS_CAi))
}






