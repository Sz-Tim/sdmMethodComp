# Auxilliary functions
# Comparison of SDM approaches
# Tim Szewczyk


##-- antilogit
antilogit <- function (x) {
  exp(x)/(1 + exp(x))
}



##-- logit
logit <- function (x) {
  log(x/(1-x))
}



##-- generate landscape
build_landscape <- function(f, x_max=Inf, y_max=Inf) {
  library(tidyverse)
  # load GIS data
  lc.df <- read_csv(f) %>% 
    filter(!is.na(bio1_mean)) %>%
    mutate(x=as.integer(factor(.$left)),
           y=as.integer(factor(.$top, levels=rev(levels(factor(.$top))))),
           x_y=paste(x, y, sep="_"))
  # establish grid with x=1:ncol, y=1:nrow
  env.rct <- as.tibble(expand.grid(x=1:max(lc.df$x), y=1:max(lc.df$y))) %>%
    mutate(x_y=paste(x, y, sep="_")) 
  match_id <- match(env.rct$x_y, lc.df$x_y)
  # pair environmental variables
  env.rct <- env.rct %>%
    mutate(temp=c(scale(lc.df$bio1_mean))[match_id],
           temp2=temp^2,
           prec=c(scale(lc.df$bio12_mean))[match_id],
           prec2=prec^2,
           pOpn=lc.df$nlcd1_mean[match_id],
           pOth=lc.df$nlcd2_mean[match_id],
           pDec=lc.df$nlcd3_mean[match_id],
           pEvg=lc.df$nlcd4_mean[match_id],
           pMxd=lc.df$nlcd5_mean[match_id],
           inbd=!is.na(match(.$x_y, lc.df$x_y)),
           lat=lc.df$top[match_id],
           lon=lc.df$left[match_id],
           rdLen=lc.df$rdLen[match_id]) %>%
    filter(x <= x_max & y >= (max(.$y)-y_max)) %>%
    mutate(id=row_number(), 
           id.inbd=min_rank(na_if(inbd*id, 0)))
  env.rct[is.na(env.rct)] <- 0
  # subset inbound cells
  env.in <- filter(env.rct, inbd) %>% select(4:12, 1:3, 13:18) %>%
    mutate(pOpn=c(scale(pOpn)), pOth=c(scale(pOth)), pDec=c(scale(pDec)),
           pEvg=c(scale(pEvg)), pMxd=c(scale(pMxd)))
  return(list(env.rct=env.rct, env.in=env.in))
}



##-- create design matrix from vector of sizes
z_pow <- function(z.vec, n_z) {
  if(n_z==1) {
    return(matrix(1, ncol=1, nrow=length(z.vec)))
  } else {
    z <- matrix(z.vec, nrow=length(z.vec), ncol=n_z-1)
    for(i in 1:(n_z-1)) z[,i] <- z[,i]^i
    z <- cbind(1, z)
    return(z)
  }
}



##-- calculate means and sds from multiple CA simulations
summarize_CA_simulations <- function(sim.ls, tmax, y.ad) {
  library(tidyverse)
  s.a <- list(B=map(sim.ls, ~.$B),
              nSd=map(sim.ls, ~.$nSd),
              nSdStay=map(sim.ls, ~.$nSdStay),
              D=map(sim.ls, ~.$D),
              N_tot=map(sim.ls, ~apply(.$N, 1:2, sum)),
              N_ad=map(sim.ls, ~.$N[,,y.ad]),
              N_rcr=map(sim.ls, ~.$N[,,1])) %>%
    map(simplify2array)
  return(list(B.mn=apply(s.a$B, 1:2, mean), 
              nSd.mn=apply(s.a$nSd, 1:2, mean), 
              nSdStay.mn=apply(s.a$nSdStay, 1:2, mean), 
              D.mn=apply(s.a$D, 1:2, mean), 
              N_tot.mn=apply(s.a$N_tot, 1:2, mean),
              N_ad.mn=apply(s.a$N_ad, 1:2, mean),
              N_rcr.mn=apply(s.a$N_rcr, 1:2, mean)))
}



##-- calculate means from multiple IPM simulations
summarize_IPM_simulations <- function(sim.ls, tmax) {
  library(tidyverse)
  s.a <- list(B=map(sim.ls, ~.$B),
              nSd=map(sim.ls, ~.$nSd),
              D=map(sim.ls, ~.$D),
              p_est=map(sim.ls, ~.$p_est.i),
              N_sim=map(sim.ls, ~.$N_sim),
              N_tot=map(sim.ls, 
                        ~map_dbl(.$d, ~sum(!is.na(.$sizeNext[.$yr==tmax])))),
              N_surv=map(sim.ls, 
                         ~map_dbl(.$d, ~sum(.$surv[.$yr==tmax], na.rm=T))),
              N_rcr=map(sim.ls, 
                        ~map_dbl(.$d, ~sum(is.na(.$size[.$yr==tmax]))))) %>% 
    map(simplify2array)
  return(list(B.mn=apply(s.a$B, 1:2, mean), 
              nSd.mn=apply(s.a$nSd, 1:2, mean), 
              D.mn=apply(s.a$D, 1:2, mean), 
              p_est.mn=apply(s.a$p_est, 1:2, mean), 
              N_sim.mn=apply(s.a$N_sim, 1:2, mean), 
              N_tot.mn=apply(s.a$N_tot, 1, mean),
              N_surv.mn=apply(s.a$N_surv, 1, mean),
              N_rcr.mn=apply(s.a$N_rcr, 1, mean)))
}



##-- calculate means from simulation means of multiple CA samples
summarize_CA_samples <- function(CA.f, in.id) {
  Sa <- list(B.mn=map(CA.f, ~.$B.mn),
             nSd.mn=map(CA.f, ~.$nSd.mn),
             nSdStay.mn=map(CA.f, ~.$nSdStay.mn),
             D.mn=map(CA.f, ~.$D.mn),
             N_tot.mn=map(CA.f, ~.$N_tot.mn),
             N_ad.mn=map(CA.f, ~.$N_ad.mn),
             N_rcr.mn=map(CA.f, ~.$N_rcr.mn)) %>%
    map(simplify2array)
  return(list(B.mn=apply(Sa$B.mn[in.id,,], 1:2, mean),
             nSd.mn=apply(Sa$nSd.mn[in.id,,], 1:2, mean),
             nSdStay.mn=apply(Sa$nSdStay.mn[in.id,,], 1:2, mean),
             D.mn=apply(Sa$D.mn[in.id,,], 1:2, mean),
             N_tot.mn=apply(Sa$N_tot.mn[in.id,,], 1:2, mean),
             N_ad.mn=apply(Sa$N_ad.mn[in.id,,], 1:2, mean),
             N_rcr.mn=apply(Sa$N_rcr.mn[in.id,,], 1:2, mean)))
}



##-- calculate means from simulation means of multiple IPM samples
summarize_IPM_samples <- function(U.f, S.f) {
  Ua <- list(IPM.mn=map(U.f, ~.$IPMs),
             P.mn=map(U.f, ~.$Ps),
             F.mn=map(U.f, ~.$Fs),
             Nt.mn=map(U.f, ~.$Nt),
             lam.mn=map(U.f, ~.$lam.t),
             p_est.mn=map(U.f, ~.$p_est.t)) %>% 
    map(simplify2array)
  Uf <- list(IPM.mn=apply(Ua$IPM.mn, 1:3, mean),
             P.mn=apply(Ua$P.mn, 1:3, mean),
             F.mn=apply(Ua$F.mn, 1:3, mean),
             Nt.mn=apply(Ua$Nt.mn, 1:3, mean),
             lam.mn=apply(Ua$lam.mn, 1:2, mean),
             p_est.mn=apply(Ua$p_est.mn, 1:2, mean))
  Sa <- list(B.mn=map(S.f, ~.$B.mn),
             nSd.mn=map(S.f, ~.$nSd.mn),
             D.mn=map(S.f, ~.$D.mn),
             p_est.mn=map(S.f, ~.$p_est.mn),
             N_sim.mn=map(S.f, ~.$N_sim.mn),
             N_tot.mn=map(S.f, ~.$N_tot.mn),
             N_surv.mn=map(S.f, ~.$N_surv.mn),
             N_rcr.mn=map(S.f, ~.$N_rcr.mn)) %>%
    map(simplify2array)
  Sf <- list(B.mn=apply(Sa$B.mn, 1:2, mean),
             nSd.mn=apply(Sa$nSd.mn, 1:2, mean),
             D.mn=apply(Sa$D.mn, 1:2, mean),
             p_est.mn=apply(Sa$p_est.mn, 1:2, mean),
             N_sim.mn=apply(Sa$N_sim.mn, 1:2, mean),
             N_tot.mn=apply(Sa$N_tot.mn, 1, mean),
             N_surv.mn=apply(Sa$N_surv.mn, 1, mean),
             N_rcr.mn=apply(Sa$N_rcr.mn, 1, mean))
  return(list(Uf=Uf, Sf=Sf))
}



##-- extract SDM sampling & modeling issues from file names
extract_SDM_issues <- function(f, SDM) {
  library(stringr)
  issues <- str_remove(str_split(f, paste0(SDM, "_"), Inf, T)[,2], ".rds")
  return(str_split(issues, "_"))
}
