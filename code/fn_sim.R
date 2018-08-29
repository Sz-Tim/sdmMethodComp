# Simulation functions
# Comparison of SDM approaches
# Tim Szewczyk



##-- simulate expected values
#' Generate expected values for simulation based on the vital rate regressions
#' from the IPM and the current size distribution
#' @param k Current year
#' @param z.ik Current size distribution in cell i
#' @param e.ik Row numbers in \code{E_i} to generate expected values for
#' @param E_i List of expected values for each individual in cell i
#' @param p List of parameters
#' @param X.i List of covariates, with elements \code{.$s, .$g, .$fl, .$seed}
#' @param n_z List of maximum exponents to raise the size distribution to
#' @return List E_i with new expected values appended
sim_expected <- function(k, z.ik, e.ik, E_i, p, X.i, n_z) {
  if(length(z.ik) > 0) {
    z.mx <- cbind(1, z.ik, z.ik^2, z.ik^3)
    E_i$yr[e.ik] <- k
    E_i$size[e.ik] <- z.ik
    E_i$s[e.ik] <- antilogit(z.mx[,1:n_z$s] %*% p$s_z + c(X.i$s %*% p$s_x))
    E_i$g[e.ik] <- z.mx[,1:n_z$g] %*% p$g_z + c(X.i$g %*% p$g_x)
    E_i$fl[e.ik] <- antilogit(z.mx[,1:n_z$fl] %*% p$fl_z + c(X.i$fl %*% p$fl_x))
    E_i$seed[e.ik] <- exp(z.mx[,1:n_z$seed] %*% p$seed_z + c(X.i$seed %*% p$seed_x))
  }
  return(E_i)
}




##-- simulate realized values
#' Generate realized values for simulation based on the expected values 
#' produced by \link{sim_expected}
#' @param k Current year
#' @param z.ik Current size distribution
#' @param d_i List of simulated data for cell i
#' @param d.ik Row numbers in \code{d_i} to generate realized values for
#' @param e.ik Row numbers in \code{E_i} to generate expected values for
#' @param E_i List of expected values for each individual in cell i
#' @param p List of parameters
#' @param lo Minimum allowable size
#' @param hi Maximum allowable size
#' @return List d_i with new realized values appended
sim_realized <- function(k, z.ik, d_i, d.ik, e.ik, E_i, p, lo, hi) {
  n_t <- length(d.ik)
  if(n_t > 0) {
    d_i$yr[d.ik] <- k
    d_i$size[d.ik] <- z.ik
    d_i$surv[d.ik] <- rbinom(n_t, 1, E_i$s[e.ik])
    surv.ik <- ifelse(d_i$surv[d.ik], 1, NA)
    d_i$sizeNext[d.ik] <- rnorm(n_t, E_i$g[e.ik], p$g_sig) * surv.ik
    d_i$sizeNext[d.ik] <- pmax(d_i$sizeNext[d.ik], lo)
    d_i$sizeNext[d.ik] <- pmin(d_i$sizeNext[d.ik], hi)
    d_i$fl[d.ik] <- fl.ik <- rbinom(n_t, 1, E_i$fl[e.ik]) * surv.ik
    d_i$seed[d.ik] <- rpois(n_t, E_i$seed[e.ik]) * ifelse(fl.ik, 1, NA)
  }
  return(d_i)
}




##-- simulate number of recruits
#' Generate size distribution for new recruits in cell i in year k, given the 
#' number of seeds produced, the size of the seed bank, and the number of 
#' immigrant seeds.
#' @param k Current year
#' @param d_i List of simulated data for cell i
#' @param p_est_ik Establishment probability for cell i in year k
#' @param nSd_ik Number of seeds produced in cell i in year k
#' @param B_ik Number of seeds in the seed bank in cell i in year k
#' @param D_ik Nmuber of immigrant seeds to cell i in year k
#' @param p List of parameters
#' @param ldd \code{FALSE} Is this to implement a long distance dispersal event?
#' @return List d_i with new recruits appended
sim_recruits <- function(k, d_i, p_est_ik, nSd_ik, B_ik, D_ik, p, ldd=F) {
  if(ldd) {
    n <- 1
  } else {
    n <- round(p_est_ik * (nSd_ik * (1-p$p_emig) * p$rcr_dir + 
                             B_ik * p$rcr_SB + 
                             D_ik * p$rcr_dir)) 
  }
  if(n > 0) {
    rcr_i <- (1:n)+length(d_i$yr)
    d_i$sizeNext[rcr_i] <- rnorm(n, p$rcr_z[1], p$rcr_z[2])
    d_i$yr[rcr_i] <- k
    d_i$size[rcr_i] <- d_i$surv[rcr_i] <- d_i$fl[rcr_i] <- d_i$seed[rcr_i] <- NA
    d_i$age[rcr_i] <- 0
  }
  return(d_i)
}




##-- simulate addition to seed bank
#' Generate new additions to the seed bank
#' @param nSd_ik Number of seeds produced in cell i in year k
#' @param B_ik Number of seeds in the seed bank in cell i in year k
#' @param D_ik Nmuber of immigrant seeds to cell i in year k
#' @param p List of parameters
#' @return Integer abundance of the seed bank in the next year
sim_seedbank <- function(nSd_ik, B_ik, D_ik, p) {
  round(p$s_SB * (B_ik * (1-p$rcr_SB) + 
                    nSd_ik * (1-p$p_emig) * (1-p$rcr_dir) + 
                    D_ik * (1-p$rcr_dir)))
}




##-- simulate full dataset
#' Simulate individual growth, reproduction, and dispersal across each cell in
#' the landscape for a set number of years. This draws on the sim_* functions
#' defined above.
#' @param n.cell Number of inbound cells in landscape
#' @param lo Minimum allowable size
#' @param hi Maximum allowable size
#' @param p List of parameters
#' @param X List of covariates, with elements \code{.$s, .$g, .$fl, .$seed}
#' @param n_z Maximum exponent to raise the size distribution to
#' @param sdd Short distance dispersal array of neighborhoods and probabilities;
#' generated with \link{gbPopMod::set_sdd_probs()$i}; perspective is the 
#' dispersal FROM each source cell i
#' @param sdd.j Short distance dispersal immigrant neighborhoods TO each cell j
#' @param N_init Vector of initial population sizes, length n.cell
#' @param verbose \code{FALSE} Show progress bar?
#' @return List with elements E = list with element for each cell with expected 
#' values for each individual in each year, d = list with element for each cell 
#' with realized values (i.e., data) for each individual in each year, nSd = 
#' matrix with number of seeds produced in each cell in each year, D = matrix
#' with number of immigrant seeds in each cell in each year, and p_est.i = 
#' density dependent establishment probabilities in each cell in each year
simulate_data <- function(n.cell, lo, hi, p, X, n_z, sdd, sdd.j, N_init, verbose=F) {
  library(tidyverse)
  i <- 1:n.cell
  if(verbose) pb <- txtProgressBar(min=1, max=p$tmax, style=3)
  
  # storage objects
  E <- d <- map(i, ~list())  # much faster than using dataframes
  B <- matrix(0, nrow=n.cell, ncol=(p$tmax+1))
  nSd <- D <- matrix(0, nrow=n.cell, ncol=p$tmax)
  p_est.i <- matrix(p$p_est, nrow=n.cell, ncol=p$tmax)
  z.k <- map(i, ~runif(N_init[.], p$z.rng[1], p$z.rng[2]))
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
    nSd[,k] <- vapply(i, function(x) sum(d[[x]]$seed[d.k[[x]]], na.rm=TRUE), 1)
    
    ## dispersal & density dependence
    sdd.D <- vapply(i, function(x) nSd[x,k]*p$p_emig*sdd[,,1,x], sdd[,,1,1])
    D[,k] <- vapply(sdd.j, function(x) sum(sdd.D[x]), 1)
    if(p$NDD) p_est.i[,k] <- pmin(p$NDD_n/(nSd[,k]+D[,k]), p$p_est)
    d <- lapply(i, function(x) sim_recruits(k, d[[x]], p_est.i[x,k], nSd[x,k], 
                                            B[x,k], D[x,k], p, F))
    ldd_k <- sample.int(n.cell, p$ldd)
    for(j in seq_along(ldd_k)) {
      d[[ldd_k[j]]] <- sim_recruits(k, d[[ldd_k[j]]], 1, 1, 0, 0, p, T)
    }
    B[,k+1] <- vapply(i, function(x) sim_seedbank(nSd[x,k], B[x,k], D[x,k],p),1)
    if(verbose) setTxtProgressBar(pb, k)
  }
  return(list(E=E, d=d, B=B, nSd=nSd, D=D, p_est.i=p_est.i))
}




