# Auxilliary functions
# Comparison of SDM approaches
# Tim Szewczyk


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
                N_tot=map(S.i, ~map_dbl(.$d, 
                                        ~sum(!is.na(.$sizeNext[.$yr==tmax])))),
                N_surv=map(S.i, ~map_dbl(.$d, 
                                         ~sum(.$surv[.$yr==tmax], na.rm=T))),
                N_rcr=map(S.i, ~map_dbl(.$d, 
                                        ~sum(is.na(.$size[.$yr==tmax]))))) %>% 
    map(simplify2array)
  rm(S.i)
  
  return(list(B.mn=apply(S.sum$B, c(1,2), mean), 
              B.sd=apply(S.sum$B, c(1,2), sd), 
              nSd.mn=apply(S.sum$nSd, c(1,2), mean), 
              nSd.sd=apply(S.sum$nSd, c(1,2), sd), 
              D.mn=apply(S.sum$D, c(1,2), mean), 
              D.sd=apply(S.sum$D, c(1,2), sd), 
              p_est.mn=apply(S.sum$p_est, c(1,2), mean), 
              p_est.sd=apply(S.sum$p_est, c(1,2), sd), 
              N_sim.mn=apply(S.sum$N_sim, c(1,2), mean), 
              N_sim.sd=apply(S.sum$N_sim, c(1,2), sd),
              N_tot.mn=apply(S.sum$N_tot, 1, mean),
              N_tot.sd=apply(S.sum$N_tot, 1, sd),
              N_surv.mn=apply(S.sum$N_surv, 1, mean),
              N_surv.sd=apply(S.sum$N_surv, 1, mean),
              N_rcr.mn=apply(S.sum$N_rcr, 1, mean),
              N_rcr.sd=apply(S.sum$N_rcr, 1, mean)))
}

