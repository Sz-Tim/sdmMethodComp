# Functions for imposing sampling & modeling issues
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







