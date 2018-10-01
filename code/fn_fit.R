# Functions for sampling, imposing issues, & fitting models
# Comparison of SDM approaches
# Tim Szewczyk



################################################################################
## Sampling functions
########

##-- sample for demographic CA model
#' Sample from S and summarise to population-level statistics for the specified
#' cells and years within each observed dataset
#' @param S True individual-level data produced in 1_simulateSpecies.R
#' @param Mech.sample List of cell indices for each set of samples
#' @param O_yr List with [["CA"]] containing the years to sample
#' @param prop.sampled Proportion of individuals (per cell) to sample
#' @return List with an element for each dataset, where each dataset is a list
#'   containing a dataframe "d" with population-level statistics for each cell
#'   and year, matrix "B" (dim=[n.cell, length(O_yr$CA)] with the size of the
#'   seedbank, and matrix "D" (dim=dim(B)) with the number of immigrant seeds to
#'   each cell
sample_for_CA <- function(S, Mech.sample, O_yr, prop.sampled) {
  O_CA <- vector("list", length(Mech.sample))
  for(s in seq_along(O_CA)) {
    CA.d <- CA.B <- CA.D <- vector("list", length(Mech.sample[[1]]))
    for(j in seq_along(Mech.sample[[s]])) {
      i <- Mech.sample[[s]][j]
      CA.d[[j]] <- data.frame(S$d[[i]]) %>% filter(yr %in% O_yr$CA)
      CA.d[[j]] <- sample_frac(CA.d[[j]], prop.sampled)
      if(!all(O_yr$CA %in% CA.d[[j]]$yr)) {
        missing.yr <- O_yr$CA[!O_yr$CA %in% CA.d[[j]]$yr]
        CA.d[[j]] <- CA.d[[j]] %>%
          add_row(yr=missing.yr, 
                  age=NA, sizeNext=NA, seed=NA, fl=NA, surv=NA, size=NA)
      }
      CA.d[[j]] <- CA.d[[j]] %>% group_by(yr) %>% 
        summarise(N=sum(!is.na(size)), 
                  s.N.0=sum(surv[age>2]==0, na.rm=TRUE),
                  s.N.1=sum(surv[age>2]==1, na.rm=TRUE),
                  s.M.0=sum(surv[age<3]==0, na.rm=TRUE),
                  s.M.1=sum(surv[age<3]==1, na.rm=TRUE),
                  f.0=sum(fl==0, na.rm=TRUE),
                  f.1=sum(fl==1, na.rm=TRUE),
                  mu=mean(seed, na.rm=TRUE) %>% round,
                  nSeed=sum(seed, na.rm=TRUE),
                  N.rcr=sum(is.na(size) & !is.na(sizeNext))) %>%
        mutate(mu=ifelse(is.nan(mu), 0, mu)) %>%
        ungroup() %>%
        mutate(lambda=N/lag(N,1)) %>%
        add_column(id.in=i) %>%
        full_join(env.in[i,], by="id.in")
      CA.B[[j]] <- S$B[i,tail(1:dim(S$B)[2], length(O_yr$CA))]
      CA.D[[j]] <- S$D[i,tail(1:dim(S$B)[2], length(O_yr$CA))]
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



##-- sample for IPM and individual-based CA model
#' Sample from S for the specified cells and years within each observed dataset
#' @param S True individual-level data produced in 1_simulateSpecies.R
#' @param Mech.sample List of cell indices for each set of samples
#' @param O_yr List with [["IPM"]] containing the years to sample
#' @param prop.sampled Proportion of individuals (per cell) to sample
#' @return List with an element for each dataset, where each dataset is a
#'   dataframe with individual-level statistics for each cell and year
sample_for_IPM <- function(S, Mech.sample, O_yr, prop.sampled) {
  O_IPM <- vector("list", length(Mech.sample))
  for(s in seq_along(O_IPM)) {
    IPM.d <- vector("list", length(Mech.sample[[1]]))
    for(j in seq_along(Mech.sample[[s]])) {
      i <- Mech.sample[[s]][j]
      IPM.d[[j]] <- data.frame(S$d[[i]]) %>% filter(yr %in% O_yr$IPM) %>%
        mutate(size2=size^2, size3=size^3) %>%
        add_column(id.in=i) %>%
        full_join(env.in[i,], by="id.in")
      IPM.d[[j]] <- sample_frac(IPM.d[[j]], prop.sampled)
    }
    O_IPM[[s]] <- do.call(rbind, IPM.d)
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
    n_obs <- nrow(O_IPM[[s]])
    O_IPM[[s]] <- O_IPM[[s]] %>%
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
  p.mod$p_emig <- pexp(0.5, p$sdd_rate, lower.tail=F)
  p.mod$n.ldd <- ldd
  return(p.mod)
}






################################################################################
## Fitting functions
########

##-- fit MaxEnt
#' Generate predicted distribution using MaxEnt
#' @param sp Virtual species
#' @param issue Issue name
#' @param samp.issue Issue with observed dataset
#' @param lam.df Dataframe with covariates, generated in 1_simulateSpecies.R
#' @param v Variables to include in model
#' @return List with dataframe [["P_MxE"]] containing a row for each cell and
#'   columns for cell information, mean(probability of presence), and
#'   sd(probability of presence); and [["diag"]] containing diagnostics
fit_MxE <- function(sp, issue, samp.issue, lam.df, v) {
  library(dismo); library(here); library(tidyverse); library(raster)
  path_iss <- paste0("out/maxent/", sp, "/", issue, "/")
  # load observations
  O_Mx <- readRDS(here("vs", sp, paste0("O_Mx_", samp.issue, ".rds")))
  X.Mx <- lam.df[,names(lam.df) %in% v]
  temp.rast <- vector("list", ncol(X.Mx))
  for(vi in 1:ncol(X.Mx)) {
    temp.rast[[vi]] <- rasterFromXYZ(cbind(lam.df[,c("lon","lat")], X.Mx[,vi]))
  }
  names(temp.rast) <- names(X.Mx)
  rast.Mx <- stack(temp.rast)

  bias_type="none"
  MxE.f <- MxE.p <- vector("list", length(O_Mx))
  fit.args <- c("responsecurves", "jackknife", "replicates=10", "plots",
                "removeduplicates", "nothreshold", "nohinge", "pictures",
                "noautofeature", "noprefixes", "writeplotdata",
                "outputformat=logistic")
  for(i in seq_along(O_Mx)) {
    if(!dir.exists(paste0(path_iss, i))) dir.create(paste0(path_iss, i), recursive=T)
    MxE.f[[i]] <- maxent(x=rast.Mx, p=as.matrix(lam.df[O_Mx[[i]], c("lon", "lat")]),
                        args=fit.args, path=paste0(path_iss, i))
    MxE.p[[i]] <- mean(predict(MxE.f[[i]], rast.Mx, args="outputformat=logistic"))
    d_j <- vector("list", 10)
    for (j in 0:9){   #loop for each replicate
      d <- read.csv(paste0(path_iss, i, "/species_", j, "_samplePredictions.csv"))
      d$run <- j
      dt <- d[d$Test.or.train=="train",]  # training data
      dtest <- d[d$Test.or.train=="test",]  # testing data
      dt <- dt[order(dt$Logistic.prediction, decreasing=T),]

      d_j[[j+1]] <- data.frame(Dataset=i,
                               Run=j,
                               thresh_Obs=dt$Logistic.prediction[nrow(dt)],
                               Bias=bias_type,
                               test_n=nrow(dtest),
                               train_n=nrow(dt))
    }
    thresh <- mean(do.call("rbind", d_j)$thresh_Obs)
    MxE.p[[i]]@data@values <- MxE.p[[i]]@data@values >= thresh
    gc()
  }
  # munge output
  names(MxE.p) <- 1:length(MxE.p)
  MxE.p.data <- map_dfr(MxE.p, ~.@data@values) %>% as.matrix
  P_MxE <- lam.df %>% 
    dplyr::select("x", "y", "x_y", "lat", "lon", "id", "id.in") %>%
    mutate(prP=apply(na.omit(MxE.p.data), 1, mean),
           prP.sd=apply(na.omit(MxE.p.data), 1, sd))
  diagnostics <- NULL
  # diagnostics <- map(MxE.f, ~evaluate(p=S_p, a=S_a, model=., x=rast.Mx))
  return(list(P_MxE=P_MxE, diag=diagnostics))
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
#' @param sp Virtual species
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
fit_CA <- function(sp, samp.issue, mod.issue, p, env.rct, env.rct.unsc, lam.df, 
                   v, m, N_init, sdd.pr, sdd.ji, p.ji, n_sim, n_cores) {
  library(here); library(tidyverse); library(gbPopMod); 
  library(MuMIn); library(lme4); library(doSNOW)
  out.dir <- paste0("out/", sp, "/", samp.issue, "/", mod.issue, "/temp/")
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  n.cell <- nrow(lam.df)
  n.grid <- nrow(env.rct)
  
  # load observations
  O_CA <- readRDS(here("vs", sp, paste0("O_CA_", samp.issue, ".rds")))
  X.CA <- env.rct %>% rename(id.in=id.in) %>%
    dplyr::select(one_of("x", "y", "x_y", "inbd", "id", "id.in", names(v)[-1]))
  n_LC <- n_distinct(read.csv(paste0("data/PNAS_2017/aggLC_", 
                                     ifelse(grepl("mustard", sp), 
                                            "mustard.csv", "woody.csv")))$agg)
  
  # initialize CA parameters
  p.CA <- set_g_p(tmax=p$tmax, 
                  lc.r=diff(range(lam.df$y)), lc.c=diff(range(lam.df$x)),
                  n.lc=5, N.p.t0=n.cell, 
                  m=3, gamma=1, 
                  s.B=p$s_SB, g.D=p$rcr_dir, g.B=p$rcr_SB,
                  p=as.matrix(logit(p$p_est)), 
                  sdd.max=p$sdd_max, sdd.rate=p$sdd_rate, n.ldd=1,
                  p.c=matrix(1), bird.hab=p$bird_hab, s.c=1, method="lm")
  p.CA$p_emig <- p$p_emig
  
  # impose dispersal issues
  if(mod.issue=="underDisp") {
    p.CA <- add_misDisperse(p.CA, p, sdd_max_adj=-2, sdd_rate_adj=2, ldd=0)
  } else if(mod.issue=="overDisp") {
    p.CA <- add_misDisperse(p.CA, p, sdd_max_adj=2, sdd_rate_adj=.5, ldd=5)
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
  foreach(i=1:length(O_CA), 
          .packages=c("tidyverse", "gbPopMod", "MuMIn", "lme4")) %dopar% {
    walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), source)
    O_CA.i <- O_CA[[i]]$d 
    O_CA.K <- O_CA.i %>% filter(id %in% O_CA.i$id[abs(O_CA.i$lambda-1)<0.05])
    sim.ls <- vector("list", n_sim)
    
    # global models
    options(na.action="na.fail")
    full.m <- opt.m <- vars.opt <- setNames(vector("list", 5), 
                                            c("K", "s.M", "s.N", "p.f", "mu"))
    
    full.m$K <- glmer(as.formula(paste("N ~", m, collapse="")),
                      data=O_CA.K, family="poisson")
    if(sum(O_CA.i$s.M.0==0)/nrow(O_CA.i) < 0.9) { # in case very low mortality
      full.m$s.M <- glmer(as.formula(paste("cbind(s.M.1, s.M.0) ~", m, collapse="")), 
                          data=O_CA.i, family="binomial")
    }
    if(sum(O_CA.i$s.N.0==0)/nrow(O_CA.i) < 0.9) { # in case very low mortality
      full.m$s.N <- glmer(as.formula(paste("cbind(s.N.1, s.N.0) ~", m, collapse="")), 
                          data=O_CA.i, family="binomial")
    } 
    full.m$p.f <- glmer(as.formula(paste("cbind(f.1, f.0) ~", m, collapse="")), 
                        data=O_CA.i, family="binomial")
    full.m$mu <- glmer(as.formula(paste("mu ~", m, collapse="")), 
                       data=O_CA.i, family="poisson")
    
    # store coefficients from optimal models
    vars.ls <- rep(list(v), length(full.m)); names(vars.ls) <- names(opt.m)
    for(j in seq_along(full.m)) {
      if(is.null(full.m[[j]])) {
        vars.opt[[j]] <- c("(Intercept)"=logit(0.9999))
      } else {
        opt.m[[j]] <- get.models(dredge(full.m[[j]], m.max=6), subset=1)[[1]]
        vars.opt[[j]] <- colMeans(coef(opt.m[[j]])$yr)
      }
      vars.ls[[j]][names(vars.opt[[j]])] <- vars.opt[[j]]
    }
    
    # update parameters
    p.CA$K <- vars.ls$K
    p.CA$s.M <- vars.ls$s.M
    p.CA$s.N <- vars.ls$s.N
    p.CA$p.f <- vars.ls$p.f
    p.CA$mu <- vars.ls$mu
    
    # impose seed bank issue
    if(mod.issue=="noSB") p.CA$s.sb <- 0
    
    # run simulations
    N.init <- matrix(0, n.grid, p.CA$m)  # column for each age class
    N.init[X.CA$id[X.CA$inbd], p.CA$m] <- N_init
    N.init[X.CA$id[X.CA$inbd], -p.CA$m] <- round(N_init/5)
    for(s in 1:n_sim) {
      sim.ls[[s]] <- gbPopMod::run_sim(n.grid, n.cell, p.CA, X.CA, sdd.pr, 
                                       N.init, NULL, F, (-1:0)+p.CA$tmax)
    }
    
    i_pad <- str_pad(i, 2, pad="0")
    saveRDS(aggregate_CAd_simulations(sim.ls, max(p.CA$m)), 
            paste0(out.dir, "/CAd_fit_", i_pad, ".rds"))
    saveRDS(calc_lambda(p.CA, X.CA, sdd.ji, p.ji, method="lm")$lambda, 
            paste0(out.dir, "/CAd_lam_", i_pad, ".rds"))
    saveRDS(list(p.CA, vars.ls), paste0(out.dir, "/CAd_diag_", i_pad, ".rds"))
  }
  stopCluster(p.c)
  
  out <- list.files(out.dir, "CAd_fit", full.names=T) %>% map(readRDS) %>%
    summarize_CAd_samples(., lam.df$id)
  lambdas <- list.files(out.dir, "CAd_lam", full.names=T) %>% map_dfc(., readRDS)
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
  P_CAl <- lam.df %>%
    dplyr::select("x", "y", "x_y", "lat", "lon", "id", "id.in") %>% 
    mutate(lambda.f=rowMeans(lambdas)[lam.df$id.in],
           prP=rowMeans(lambdas>1)[lam.df$id.in])
  
  return(list(P_CAd=P_CAd,  P_CAl=P_CAl, diag=diagnostics))
}



##-- fit IPM & individual-based CA model
#' Generate predicted distribution using an IPM and an individual-based CA model
#' @param sp Virtual species
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
fit_IPM <- function(sp, samp.issue, mod.issue, p, env.rct.unsc, lam.df, v, m, 
                    n_z, n_x, N_init, sdd.ji, p.ji, n_sim, n_cores) {
  library(here); library(tidyverse); library(magrittr); library(gbPopMod);
  library(MuMIn); library(doSNOW)
  walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), source)
  out.dir <- paste0("out/", sp, "/", samp.issue, "/", mod.issue, "/temp/")
  if(!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
  n.cell <- nrow(lam.df)
  
  # load observations
  O_IPM <- readRDS(here("vs", sp, paste0("O_IPM_", samp.issue, ".rds")))
  X.IPM <- map(n_x, ~as.matrix(lam.df[,names(lam.df) %in% names(v)[-(1:n_z$s)]]))
  
  # initialize IPM parameters
  p.IPM <- p
  
  # impose dispersal issues
  if(mod.issue=="underDisp") {
    p.IPM <- add_misDisperse(p.IPM, p, sdd_max_adj=-2, sdd_rate_adj=2, ldd=0)
  } else if(mod.issue=="overDisp") {
    p.IPM <- add_misDisperse(p.IPM, p, sdd_max_adj=2, sdd_rate_adj=.5, ldd=5)
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
  foreach(i=1:length(O_IPM), 
          .packages=c("tidyverse", "magrittr", "gbPopMod", "MuMIn")) %dopar% {
    walk(paste0("code/fn_", c("aux", "sim", "IPM", "fit"), ".R"), source)
    sim.ls <- vector("list", n_sim)
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
    opt.m <- map(full.m, ~get.models(dredge(., m.lim=c(1,6)), subset=1)[[1]])
    vars.opt <- map(opt.m, coef)
    vars.ls <- rep(list(v), 4); names(vars.ls) <- names(opt.m)
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
    p.IPM$rcr_z <- filter(O_IPM[[i]], is.na(size)) %>% 
      summarise(mn=mean(sizeNext), sd=sd(sizeNext)) %>% unlist
    
    # impose issues
    if(mod.issue=="noSB") p.IPM$s_SB <- 0
    
    # use estimated slopes to fill IPM matrix
    U.f <- fill_IPM_matrices(n.cell, buffer=0, discrete=1, p.IPM, n_z,
                                  n_x, X.IPM, sdd.ji, p.ji)
    
    # use estimated slopes to generate simulated data
    for(s in 1:n_sim) {
      sim.ls[[s]] <- simulate_data(n.cell, U.f$lo, U.f$hi, p.IPM, X.IPM, n_z, 
                                   sdd.ji, p.ji, N_init, save_yrs=p.IPM$tmax)
    }

    i_pad <- str_pad(i, 2, pad="0")
    saveRDS(aggregate_CAi_simulations(sim.ls, p.IPM$tmax), 
            paste0(out.dir, "/CAi_fit_", i_pad, ".rds"))
    saveRDS(U.f, paste0(out.dir, "/IPM_fit_", i_pad, ".rds"))
    saveRDS(list(p.IPM, vars.ls), paste0(out.dir, "/IPM_diag_", i_pad, ".rds"))
  }
  stopCluster(p.c)
  
  out <- summarize_IPM_CAi_samples(
    map(list.files(out.dir, "IPM_fit", full.names=T), readRDS),
    map(list.files(out.dir, "CAi_fit", full.names=T), readRDS)
  )
  diagnostics <- list.files(out.dir, "CAd_diag", full.names=T) %>% map(readRDS)
  
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
  
  return(list(P_IPM=P_IPM, P_CAi=P_CAi, diag=diagnostics))
}






