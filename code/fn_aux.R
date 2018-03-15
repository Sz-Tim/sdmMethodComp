# Auxilliary functions
# Comparison of SDM approaches
# Tim Szewczyk


##-- generate landscape
build_landscape <- function(f, x_max=Inf, y_max=Inf) {
  require(tidyverse)
  
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


##-- calculate means and sds from multiple IPM simulations
summarize_IPM_simulations <- function(S.i) {
  
  require(tidyverse)
  
  S.sum <- list(B=map(S.i, ~.$B),
                nSd=map(S.i, ~.$nSd),
                D=map(S.i, ~.$D),
                p_est=map(S.i, ~.$p_est.i),
                N_sim=map(S.i, ~.$N_sim),
                N_tot=map(S.i, 
                          ~map_dbl(.$d, ~sum(!is.na(.$sizeNext[.$yr==tmax])))),
                N_surv=map(S.i, 
                           ~map_dbl(.$d, ~sum(.$surv[.$yr==tmax], na.rm=T))),
                N_rcr=map(S.i, 
                          ~map_dbl(.$d, ~sum(is.na(.$size[.$yr==tmax]))))) %>% 
    map(simplify2array)
  rm(S.i)
  
  return(list(B.mn=apply(S.sum$B, 1:2, mean), 
              B.sd=apply(S.sum$B, 1:2, sd), 
              nSd.mn=apply(S.sum$nSd, 1:2, mean), 
              nSd.sd=apply(S.sum$nSd, 1:2, sd), 
              D.mn=apply(S.sum$D, 1:2, mean), 
              D.sd=apply(S.sum$D, 1:2, sd), 
              p_est.mn=apply(S.sum$p_est, 1:2, mean), 
              p_est.sd=apply(S.sum$p_est, 1:2, sd), 
              N_sim.mn=apply(S.sum$N_sim, 1:2, mean), 
              N_sim.sd=apply(S.sum$N_sim, 1:2, sd),
              N_tot.mn=apply(S.sum$N_tot, 1, mean),
              N_tot.sd=apply(S.sum$N_tot, 1, sd),
              N_surv.mn=apply(S.sum$N_surv, 1, mean),
              N_surv.sd=apply(S.sum$N_surv, 1, mean),
              N_rcr.mn=apply(S.sum$N_rcr, 1, mean),
              N_rcr.sd=apply(S.sum$N_rcr, 1, mean)))
}


##-- antilogit
antilogit <- function (x) {
  exp(x)/(1 + exp(x))
}
