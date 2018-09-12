# Simulation functions
# Comparison of SDM approaches
# Tim Szewczyk



##-- simulate expected values
#' Generate expected values for simulation based on the vital rate regressions
#' from the IPM and the current size distribution
#' @param k Current year
#' @param z.ik Current size distribution in cell i
#' @param p List of parameters
#' @param X.i List of covariates, with elements \code{.$s, .$g, .$fl, .$seed}
#' @param n_z List of maximum exponents to raise the size distribution to
#' @return List E_ik of expected values
sim_expected <- function(k, z.ik, p, X.i, n_z) {
  E_ik <- setNames(map(1:6, ~numeric(0L)), 
                   c("yr", "size", "s", "g", "fl", "seed"))
  n.ik <- length(z.ik)
  if(n.ik > 0) {
    z.mx <- cbind(1, z.ik, z.ik^2, z.ik^3)
    E_ik$yr <- rep(k, n.ik)
    E_ik$size <- z.ik
    E_ik$s <- antilogit(z.mx[,1:n_z$s] %*% p$s_z + c(X.i$s %*% p$s_x))
    E_ik$g <- z.mx[,1:n_z$g] %*% p$g_z + c(X.i$g %*% p$g_x)
    E_ik$fl <- antilogit(z.mx[,1:n_z$fl] %*% p$fl_z + c(X.i$fl %*% p$fl_x))
    E_ik$seed <- exp(z.mx[,1:n_z$seed] %*% p$seed_z + c(X.i$seed %*% p$seed_x))
  }
  return(E_ik)
}




##-- simulate realized values
#' Generate realized values for simulation based on the expected values 
#' produced by \link{sim_expected}
#' @param k Current year
#' @param z.ik Current size distribution
#' @param E_ik List of expected values for each individual in cell i in year k
#' @param p List of parameters
#' @param lo Minimum allowable size
#' @param hi Maximum allowable size
#' @return List d_ik with new realized values appended
sim_realized <- function(k, z.ik, E_ik, p, lo, hi) {
  d_ik <- setNames(map(1:7, ~numeric(0L)), 
                   c("yr", "size", "surv", "sizeNext", "fl", "seed", "age"))
  n.ik <- length(z.ik)
  if(n.ik > 0) {
    d_ik$yr <- rep(k, n.ik)
    d_ik$size <- z.ik
    d_ik$surv <- rbinom(n.ik, 1, E_ik$s)
    surv.ik <- ifelse(d_ik$surv, 1, NA)
    d_ik$sizeNext <- rnorm(n.ik, E_ik$g, p$g_sig) * surv.ik
    d_ik$sizeNext <- pmax(d_ik$sizeNext, lo)
    d_ik$sizeNext <- pmin(d_ik$sizeNext, hi)
    d_ik$fl <- fl.ik <- rbinom(n.ik, 1, E_ik$fl) * surv.ik
    d_ik$seed <- rpois(n.ik, E_ik$seed) * ifelse(fl.ik, 1, NA)
  }
  return(d_ik)
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
#' @param pr_germ Germination probability, if environmentally variable; if not
#'   \code{NULL}, it replaces both \code{p$rcr_SB} and \code{p$rcr_dir}
#' @param ldd \code{FALSE} Is this to implement a long distance dispersal event?
#' @return List d_i with new recruits appended
sim_recruits <- function(k, d_i, p_est_ik, nSd_ik, B_ik, D_ik, p, pr_germ, ldd=F) {
  if(!is.null(pr_germ)) p$rcr_SB <- p$rcr_dir <- pr_germ
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
#' @param pr_germ Germination probability, if environmentally variable; if not
#'   \code{NULL}, it replaces both \code{p$rcr_SB} and \code{p$rcr_dir}
#' @return Integer abundance of the seed bank in the next year
sim_seedbank <- function(nSd_ik, B_ik, D_ik, p, pr_germ) {
  if(!is.null(pr_germ)) p$rcr_SB <- p$rcr_dir <- pr_germ
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
#' @param sdd.ji Dispersing cells TO each target cell j
#' @param p.ji Dispersal probabilities TO each target cell j
#' @param N_init Vector of initial population sizes, length n.cell
#' @param save_yrs \code{"all"} Vector of years to store, or "all" to store all
#' @param verbose \code{FALSE} Show progress bar?
#' @return List with elements E = list with element for each cell with expected 
#' values for each individual in each year, d = list with element for each cell 
#' with realized values (i.e., data) for each individual in each year, nSd = 
#' matrix with number of seeds produced in each cell in each year, D = matrix
#' with number of immigrant seeds in each cell in each year, and p_est.i = 
#' density dependent establishment probabilities in each cell in each year
simulate_data <- function(n.cell, lo, hi, p, X, n_z, sdd.ji, p.ji, N_init, 
                          save_yrs="all", verbose=F) {
  library(tidyverse)
  i <- 1:n.cell
  
  # identify which years to store data
  if(save_yrs[1]=="all") {
    yrs <- setNames(1:p$tmax, 1:p$tmax)
  } else {
    yrs <- setNames(1:length(save_yrs), save_yrs)
  }
  
  # does germination probability vary with environment?
  if(!is.null(X$germ)) { 
    pr_germ <- c(antilogit(X$germ %*% p$germ_x))
  } else { 
    pr_germ <- NULL 
  }
  
  # storage objects
  d <- map(i, ~map(1:7, ~numeric(0L)) %>% 
             setNames(c("yr", "size", "surv", "sizeNext", "fl", "seed", "age")))
  E <- map(i, ~map(1:6, ~numeric(0L)) %>% 
             setNames(c("yr", "size", "s", "g", "fl", "seed")))
  B <- nSd <- D <- matrix(0, nrow=n.cell, ncol=length(yrs))
  p_est.i <- matrix(p$p_est, nrow=n.cell, ncol=length(yrs))
  
  # initial size distribution, seed bank, & covariates
  z.k <- map(i, ~runif(N_init[.], p$z.rng[1], p$z.rng[2]))
  B1 <- rep(0, n.cell)
  X_map <- lapply(i, function(x) map(X, ~.[x,]))
  
  if(verbose) pb <- txtProgressBar(min=1, max=p$tmax, style=3)
  for(k in 1:p$tmax) {
    ## local growth
    if(k>1) {
      z.k <- map(d1, ~.$sizeNext[!is.na(.$sizeNext)])
      age.k <- map(d1, ~.$age[!is.na(.$sizeNext)] + 1)
    }
    E1 <- lapply(i, function(x) sim_expected(k, z.k[[x]], p, X_map[[x]], n_z))
    d1 <- lapply(i, function(x) sim_realized(k, z.k[[x]], E1[[x]], p, lo, hi))
    if(k==1) { 
      invisible(lapply(i, function(x) d1[[x]]$age <<- rep(1, length(z.k[[x]]))))
    } else { 
      invisible(lapply(i, function(x) d1[[x]]$age <<- age.k[[x]])) 
    }
    nSd1 <- vapply(i, function(x) sum(d1[[x]]$seed, na.rm=TRUE), 1)
    
    ## dispersal & density dependence
    D1 <- sapply(i, function(x) sum(nSd1[sdd.ji[[x]]] * p$p_emig * p.ji[[x]]))
    if(p$NDD) p_est1 <- pmin(p$NDD_n/(nSd1+D1), p$p_est)
    d1 <- lapply(i, function(x) sim_recruits(k, d1[[x]], p_est1[x], nSd1[x], 
                                            B1[x], D1[x], p, pr_germ[x], F))
    ldd_k <- sample.int(n.cell, p$ldd)
    for(j in seq_along(ldd_k)) {
      d1[[ldd_k[j]]] <- sim_recruits(k, d1[[ldd_k[j]]], 1, 1, 0, 0, p, NULL, T)
    }
    B1 <- vapply(i, function(x) sim_seedbank(nSd1[x], B1[x], D1[x], 
                                             p, pr_germ[x]), 1)
    
    # store as necessary
    if(any(k %in% names(yrs))) {
      yr_k <- which(names(yrs)==k)
      d <- map2(d, d1, ~mapply(append, .x, .y, SIMPLIFY=FALSE))
      E <- map2(E, E1, ~mapply(append, .x, .y, SIMPLIFY=FALSE))
      B[,yr_k] <- B1
      nSd[,yr_k] <- nSd1
      D[,yr_k] <- D1
      p_est.i[,yr_k] <- p_est1
    }
    if(verbose) setTxtProgressBar(pb, k)
  }
  return(list(E=E, d=d, B=B, nSd=nSd, D=D, p_est.i=p_est.i))
}




