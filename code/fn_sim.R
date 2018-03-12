



##--- simulate expected values
## 
sim_expected <- function(d_i, d.k, e.k, E, p, X, n_z) {
  size.mx <- with(d_i[d.k,], cbind(1, size, size^2, size^3))
  E$s[e.k] <- antilogit(size.mx[,1:n_z$s] %*% p$s_z + c(X$s %*% p$s_x))
  E$g[e.k] <- size.mx[,1:n_z$g] %*% p$g_z + c(X$g %*% p$g_x)
  E$fl[e.k] <- antilogit(size.mx[,1:n_z$fl] %*% p$fl_z + c(X$fl %*% p$fl_x))
  E$seed[e.k] <- exp(size.mx[,1:n_z$seed] %*% p$seed_z + c(X$seed %*% p$seed_x))
  return(E)
}


##--- simulate realized values
##
sim_realized <- function(d_i, d.k, e.k, E, p, lo, hi) {
  n_t <- length(d.k)
  d_i$surv[d.k] <- rbinom(n_t, 1, E$s[e.k])
  d_i$sizeNext[d.k] <- rnorm(n_t, E$g[e.k], p$g_sig) * 
    ifelse(d_i$surv[d.k], 1, NA)
  d_i$sizeNext[d.k] <- ifelse(d_i$sizeNext[d.k] < lo, 
                                  lo, d_i$sizeNext[d.k])
  d_i$sizeNext[d.k] <- ifelse(d_i$sizeNext[d.k] > hi, 
                                  hi, d_i$sizeNext[d.k])
  d_i$fl[d.k] <- rbinom(n_t, 1, E$fl[e.k]) * ifelse(d_i$surv[d.k], 1, NA)
  d_i$seed[d.k] <- rpois(n_t, E$seed[e.k]) * 
    ifelse(d_i$surv[d.k] & d_i$fl[d.k], 1, NA)
  return(d_i)
}


##--- simulate number of recruits
##
sim_recruits <- function(d_i, p_est_ik, nSd_ik, B_ik, D_ik, p) {
  n <- round(p_est_ik * (nSd_ik * (1-p$p_emig) * p$rcr_dir + 
                           B_ik * p$rcr_SB + 
                           D_ik * p$rcr_dir)) 
  if(n > 0) {
    n_index <- (1:n)+nrow(d_i)
    d_i[n_index, "sizeNext"] <- rnorm(n, p$rcr_z[1], p$rcr_z[2])
    d_i[n_index, "yr"] <- max(d_i$yr, na.rm=T)
  }
  return(d_i)
}


##--- simulate addition to seed bank
##
sim_seedbank <- function(nSd_ik, B_ik, D_ik, p) {
  p$s_SB * (B_ik * (1-p$rcr_SB) + 
              nSd_ik * (1-p$p_emig) * (1-p$rcr_dir) + 
              D_ik * (1-p$rcr_dir))
}


##--- simulate full dataset
##
simulate_data <- function(n.cell, tmax, n_i, z.rng, lo, hi, 
                          p, X, n_z, sdd, sdd.j, verbose=FALSE) {
  require(tidyverse)
  
  # storage objects
  E <- d <- map(1:n.cell, ~data.frame(NULL))
  B <- cbind(rpois(n.cell, n_i), matrix(0, nrow=n.cell, ncol=tmax))
  nSd <- D <- p_est.i <- N_sim <- matrix(0, nrow=n.cell, ncol=tmax)
  
  for(k in 1:tmax) {
    ## local growth
    for(i in 1:n.cell) {
      if(k==1) {
        size.k <- runif(n=n_i[i], min=z.rng[1], max=z.rng[2])
      } else {
        size.k <- d[[i]]$sizeNext[!is.na(d[[i]]$sizeNext) & d[[i]]$yr == k-1]
      }
      d.k <- (1:length(size.k)) + nrow(d[[i]])
      e.k <- (1:length(size.k)) + nrow(E[[i]])
      if(length(size.k) > 0) {
        E[[i]][e.k, "yr"] <- d[[i]][d.k, "yr"] <- k
        E[[i]][e.k, "size"] <- d[[i]][d.k, "size"] <- size.k
        E[[i]] <- sim_expected(d[[i]], d.k, e.k, E[[i]], p, map(X, ~.[i,]), n_z)
        d[[i]] <- sim_realized(d[[i]], d.k, e.k, E[[i]], p, lo, hi)
        nSd[i,k] <- sum(d[[i]]$seed[d.k], na.rm=T)
        N_sim[i,k] <- length(size.k)
      }
      sdd[,,3,i] <- nSd[i,k] * p$p_emig * sdd[,,1,i]
    }
    ## dispersal & density dependence
    for(i in 1:n.cell) {
      D[i,k] <- sum(sdd[,,3,][sdd.j[[i]]])
      p_est.i[i,k] <- ifelse(p$NDD, min(p$NDD_n/(nSd[i,k]+D[i,k]), p$p_est), 
                             p$p_est)
      d[[i]] <- sim_recruits(d[[i]], p_est.i[i,k], nSd[i,k], B[i,k], D[i,k], p)
      B[i,k+1] <- sim_seedbank(nSd[i,k], B[i,k], D[i,k], p)
    }
    if(verbose) cat("Finished", n.cell, "cells for year", k, "\n")
  }
  
  return(list(E=E, d=d, B=B, nSd=nSd, D=D, p_est.i=p_est.i, N_sim=N_sim))
}





