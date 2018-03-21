



##--- simulate expected values
## 
sim_expected <- function(k, z.i, e.k, E, p, X, n_z) {
  if(length(z.i) > 0) {
    z.mx <- cbind(1, z.i, z.i^2, z.i^3)
    E$yr[e.k] <- k
    E$size[e.k] <- z.i
    E$s[e.k] <- antilogit(z.mx[,1:n_z$s] %*% p$s_z + c(X$s %*% p$s_x))
    E$g[e.k] <- z.mx[,1:n_z$g] %*% p$g_z + c(X$g %*% p$g_x)
    E$fl[e.k] <- antilogit(z.mx[,1:n_z$fl] %*% p$fl_z + c(X$fl %*% p$fl_x))
    E$seed[e.k] <- exp(z.mx[,1:n_z$seed] %*% p$seed_z + c(X$seed %*% p$seed_x))
  }
  return(E)
}


##--- simulate realized values
##
sim_realized <- function(k, z.i, d_i, d.k, e.k, E, p, lo, hi) {
  n_t <- length(d.k)
  if(n_t > 0) {
    d_i$yr[d.k] <- k
    d_i$size[d.k] <- z.i
    d_i$surv[d.k] <- rbinom(n_t, 1, E$s[e.k])
    surv.ik <- ifelse(d_i$surv[d.k], 1, NA)
    d_i$sizeNext[d.k] <- rnorm(n_t, E$g[e.k], p$g_sig) * surv.ik
    d_i$sizeNext[d.k] <- pmax(d_i$sizeNext[d.k], lo)
    d_i$sizeNext[d.k] <- pmin(d_i$sizeNext[d.k], hi)
    d_i$fl[d.k] <- fl.ik <- rbinom(n_t, 1, E$fl[e.k]) * surv.ik
    d_i$seed[d.k] <- rpois(n_t, E$seed[e.k]) * ifelse(fl.ik, 1, NA)
  }
  return(d_i)
}


##--- simulate number of recruits
##
sim_recruits <- function(k, d_i, p_est_ik, nSd_ik, B_ik, D_ik, p) {
  n <- round(p_est_ik * (nSd_ik * (1-p$p_emig) * p$rcr_dir + 
                           B_ik * p$rcr_SB + 
                           D_ik * p$rcr_dir)) 
  if(n > 0) {
    rcr_i <- (1:n)+length(d_i$yr)
    d_i$sizeNext[rcr_i] <- rnorm(n, p$rcr_z[1], p$rcr_z[2])
    d_i$yr[rcr_i] <- k
    d_i$size[rcr_i] <- d_i$surv[rcr_i] <- d_i$fl[rcr_i] <- d_i$seed[rcr_i] <- NA
    d_i$age[rcr_i] <- 0
  }
  return(d_i)
}


##--- simulate addition to seed bank
##
sim_seedbank <- function(nSd_ik, B_ik, D_ik, p) {
  round(p$s_SB * (B_ik * (1-p$rcr_SB) + 
                    nSd_ik * (1-p$p_emig) * (1-p$rcr_dir) + 
                    D_ik * (1-p$rcr_dir)))
}


##--- simulate full dataset
##
simulate_data <- function(n.cell, lo, hi, p, X, n_z, sdd, sdd.j, verbose=F) {
  require(tidyverse)
  i <- 1:n.cell
  
  # storage objects
  E <- d <- map(i, ~list())  # much faster than using dataframes
  B <- cbind(rpois(n.cell, p$n0), matrix(0, nrow=n.cell, ncol=p$tmax))
  nSd <- D <- N_sim <- matrix(0, nrow=n.cell, ncol=p$tmax)
  p_est.i <- matrix(p$p_est, nrow=n.cell, ncol=p$tmax)
  z.k <- map(i, ~runif(p$n0, p$z.rng[1], p$z.rng[2]))
  X_map <- lapply(i, function(x) map(X, ~.[x,]))
  
  for(k in 1:p$tmax) {
    ## local growth
    if(k>1) {
      z.k <- map(d, ~.$sizeNext[!is.na(.$sizeNext) & .$yr == k-1])
      age.k <- map(d, ~.$age[!is.na(.$sizeNext) & .$yr == k-1] + 1)
    }
    d.k <- lapply(i, function(x) seq_along(z.k[[x]]) + length(d[[x]]$yr))
    e.k <- lapply(i, function(x) seq_along(z.k[[x]]) + length(E[[x]]$yr))
    E <- lapply(i, function(x) sim_expected(k, z.k[[x]], e.k[[x]], E[[x]], 
                                            p, X_map[[x]], n_z))
    d <- lapply(i, function(x) sim_realized(k, z.k[[x]], d[[x]], d.k[[x]], 
                                            e.k[[x]], E[[x]], p, lo, hi))
    if(k==1) { 
      invisible(lapply(i, function(x) d[[x]]$age[d.k[[x]]] <<- 1))
    } else { 
      invisible(lapply(i, function(x) d[[x]]$age[d.k[[x]]] <<- age.k[[x]])) 
    }
    N_sim[,k] <- vapply(z.k, length, 1)
    nSd[,k] <- vapply(i, function(x) sum(d[[x]]$seed[d.k[[x]]], na.rm=TRUE), 1)
    
    ## dispersal & density dependence
    sdd.D <- vapply(i, function(x) nSd[x,k]*p$p_emig*sdd[,,1,x], sdd[,,1,1])
    D[,k] <- vapply(sdd.j, function(x) sum(sdd.D[x]), 1)
    if(p$NDD) p_est.i[,k] <- pmin(p$NDD_n/(nSd[,k]+D[,k]), p$p_est)
    d <- lapply(i, function(x) sim_recruits(k, d[[x]], p_est.i[x,k], nSd[x,k], 
                                            B[x,k], D[x,k], p))
    B[,k+1] <- vapply(i, function(x) sim_seedbank(nSd[x,k], B[x,k], D[x,k],p),1)
    if(verbose) cat("Finished", n.cell, "cells for year", k, "\n")
  }
  return(list(E=E, d=d, B=B, nSd=nSd, D=D, p_est.i=p_est.i, N_sim=N_sim))
}




