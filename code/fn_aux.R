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
  # store scaling mean and sd
  X.col <- c(16, 18, 4:8)
  lc.scale <- select(lc.df, X.col) %>% scale
  scale.i <- cbind(sc_mn=attributes(lc.scale)$`scaled:center`,
                   sc_sd=attributes(lc.scale)$`scaled:scale`)
  rownames(scale.i) <- c("temp", "prec", "pOpn", "pOth", "pDec", "pEvg", "pMxd")
  # establish grid with x=1:ncol, y=1:nrow
  env.rct <- as.tibble(expand.grid(x=1:max(lc.df$x), y=1:max(lc.df$y))) %>%
    mutate(x_y=paste(x, y, sep="_")) 
  match_id <- match(env.rct$x_y, lc.df$x_y)
  # pair environmental variables
  env.rct <- env.rct %>%
    mutate(temp=lc.scale[match_id,1],
           temp2=temp^2,
           prec=lc.scale[match_id,2],
           prec2=prec^2,
           pOpn=lc.scale[match_id,3],
           pOth=lc.scale[match_id,4],
           pDec=lc.scale[match_id,5],
           pEvg=lc.scale[match_id,6],
           pMxd=lc.scale[match_id,7],
           inbd=!is.na(match(.$x_y, lc.df$x_y)),
           lat=lc.df$top[match_id],
           lon=lc.df$left[match_id],
           rdLen=lc.df$rdLen[match_id]) %>%
    filter(x <= x_max & y <= y_max) %>%
    mutate(id=row_number(), 
           id.inbd=min_rank(na_if(inbd*id, 0)))
  env.rct[is.na(env.rct)] <- 0
  match_id <- match(env.rct$x_y, lc.df$x_y)
  env.rct.unscaled <- env.rct %>%
    mutate(pOpn=lc.df$nlcd1_mean[match_id],
           pOth=lc.df$nlcd2_mean[match_id],
           pDec=lc.df$nlcd3_mean[match_id],
           pEvg=lc.df$nlcd4_mean[match_id],
           pMxd=lc.df$nlcd5_mean[match_id])
  # subset inbound cells
  env.in <- filter(env.rct, inbd) %>% select(4:12, 1:3, 13:18) 
  return(list(env.rct=env.rct, env.rct.unscaled=env.rct.unscaled,
              env.in=env.in, scale.i=scale.i))
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
aggregate_CA_simulations <- function(sim.ls, tmax, y.ad) {
  library(tidyverse)
  l <- list(P=map(sim.ls, ~.$N[,,y.ad]>0),
            B=map(sim.ls, ~.$B),
            nSd=map(sim.ls, ~.$nSd),
            nSdStay=map(sim.ls, ~.$nSdStay),
            D=map(sim.ls, ~.$D),
            N_tot=map(sim.ls, ~apply(.$N, 1:2, sum)),
            N_ad=map(sim.ls, ~.$N[,,y.ad]),
            N_rcr=map(sim.ls, ~.$N[,,1]))
  return(map(l, simplify2array))
}



##-- calculate means from multiple IPM simulations
aggregate_IPM_simulations <- function(sim.ls, tmax) {
  library(tidyverse)
  l <- list(P=map(sim.ls, ~map_dbl(.$d, ~sum(.$surv[.$yr==tmax], na.rm=T))>0),
            B=map(sim.ls, ~.$B),
            nSd=map(sim.ls, ~.$nSd),
            D=map(sim.ls, ~.$D),
            p_est=map(sim.ls, ~.$p_est.i),
            N_sim=map(sim.ls, ~.$N_sim),
            N_tot=map(sim.ls, 
                      ~map_dbl(.$d, ~sum(!is.na(.$sizeNext[.$yr==tmax])))),
            N_surv=map(sim.ls, 
                       ~map_dbl(.$d, ~sum(.$surv[.$yr==tmax], na.rm=T))),
            N_rcr=map(sim.ls, ~map_dbl(.$d, ~sum(is.na(.$size[.$yr==tmax])))))
  return(map(l, simplify2array))
}



##-- calculate means from simulation means of multiple CA samples
summarize_CA_samples <- function(CA.f, in.id) {
  Sa <- list(P=map(CA.f, ~.$P),
             B=map(CA.f, ~.$B),
             nSd=map(CA.f, ~.$nSd),
             nSdStay=map(CA.f, ~.$nSdStay),
             D=map(CA.f, ~.$D),
             N_tot=map(CA.f, ~.$N_tot),
             N_ad=map(CA.f, ~.$N_ad),
             N_rcr=map(CA.f, ~.$N_rcr)) %>%
    map(simplify2array)
  return(list(prP=apply(Sa$P[in.id,,,], 1:2, mean, na.rm=T),
              B.mn=apply(Sa$B[in.id,,,], 1:2, mean, na.rm=T),
              nSd.mn=apply(Sa$nSd[in.id,,,], 1:2, mean, na.rm=T),
              nSdStay.mn=apply(Sa$nSdStay[in.id,,,], 1:2, mean, na.rm=T),
              D.mn=apply(Sa$D[in.id,,,], 1:2, mean, na.rm=T),
              N_tot.mn=apply(Sa$N_tot[in.id,,,], 1:2, mean, na.rm=T),
              N_ad.mn=apply(Sa$N_ad[in.id,,,], 1:2, mean, na.rm=T),
              N_rcr.mn=apply(Sa$N_rcr[in.id,,,], 1:2, mean, na.rm=T)))
}



##-- calculate means from simulation means of multiple IPM samples
summarize_IPM_samples <- function(U.f, S.f) {
  Ua <- list(IPMs=map(U.f, ~.$IPMs),
             Ps=map(U.f, ~.$Ps),
             Fs=map(U.f, ~.$Fs)) %>% 
    map(simplify2array)
  Ua_lam <- apply(Ua$IPMs, 3:4, function(x) Re(eigen(x)$values[1]))
  Uf <- list(prP=apply(Ua_lam>=1, 1, mean),
             lam.mn=apply(Ua_lam, 1, mean),
             IPM.mn=apply(Ua$IPMs, 1:3, mean),
             P.mn=apply(Ua$Ps, 1:3, mean),
             F.mn=apply(Ua$Fs, 1:3, mean))
  Sa <- list(P=map(S.f, ~.$P),
             B=map(S.f, ~.$B),
             nSd=map(S.f, ~.$nSd),
             D=map(S.f, ~.$D),
             p_est=map(S.f, ~.$p_est),
             N_sim=map(S.f, ~.$N_sim),
             N_tot=map(S.f, ~.$N_tot),
             N_surv=map(S.f, ~.$N_surv),
             N_rcr=map(S.f, ~.$N_rcr)) %>%
    map(simplify2array)
  Sf <- list(prP=apply(Sa$P, 1, mean),
             B.mn=apply(Sa$B, 1:2, mean),
             nSd.mn=apply(Sa$nSd, 1:2, mean),
             D.mn=apply(Sa$D, 1:2, mean),
             p_est.mn=apply(Sa$p_est, 1:2, mean),
             N_sim.mn=apply(Sa$N_sim, 1:2, mean),
             N_tot.mn=apply(Sa$N_tot, 1, mean),
             N_surv.mn=apply(Sa$N_surv, 1, mean),
             N_rcr.mn=apply(Sa$N_rcr, 1, mean))
  return(list(Uf=Uf, Sf=Sf))
}



##-- extract SDM sampling & modeling issues from file names
extract_SDM_details <- function(f) {
  library(stringr)
  str_split(str_remove(str_split(f, "P_", Inf, T)[,2], ".rds"), "_")
}
