# IPM vital rate functions
# Comparison of SDM approaches
# Tim Szewczyk

# This script contains functions used for constructing an IPM matrix. All 
# coefficients and probabilities must be estimated separately and used as 
# inputs for these functions. These should be contained within the list `p`,
# with a scalar or vector for each set of covariates, named as:
#
# p=list(s_z=c(-1.75, .75),  # b1 + b2*z + ...
#        s_x=c(1, -.1, -.5, -.1),  # b1*x1 + ...
#        g_z=c(.2, 2, -0.1),  # b1 + b2*z + b3*z^2 + ...
#        g_x=c(1, -.1, 1.5, -.1),  #b1*x1 + ...
#        g_sig=1,  # growth ~ N(E, g_sig)
#        fl_z=c(-1.5, .1, .1),  # b1 + b2*z + b3*z^2 + b4*z^3 + ...
#        fl_x=c(1, -.01, -1, -.01),  # b1*x1 + ...
#        seed_z=c(2, 0.5, -.03),  # b1 + b2*z + b2*z^2 + ...
#        seed_x=c(1, -.1, 1, -.1),  # b1*x1 + ...
#        rcr_z=c(1.5, 0.4),  # N(mean=rcrt1, sd=rcrt2)
#        p_est=0.02,  # p(establishment)
#        rcr_SB=0.2,  # p(recruit from seedbank)
#        rcr_dir=0.6,  # p(recruit directly)
#        s_SB=0.3,  # p(survive in seedbank additional year)
#        p_emig=0.3,  # p(seed emigrates)
#        sdd_max=2,  # max SDD distance in cells
#        sdd_rate=0.8  # SDD dispersal rate
#        )
#
# The function `z_pow()` takes a vector of sizes (z) and the number of size 
# slopes for each function, returning a design matrix with the correct number 
# of columns where: z.mx[,1]=1; z.mx[,2]=z; z.mx[,3]=z^2; z.mx[,4]=z^3; etc


################################################################################
## Vital rate functions
########

##-- Survival
##   P(s) ~ antilogit(size + environment)
#' Calculate individual survival probabilities based on size & environment
#' @param z.v Vector of sizes: discretized size range in approximated IPM matrix
#' @param p List of parameters
#' @param n_sz Maximum exponent to raise the size distribution to
#' @param X.s Matrix of environmental covariates 
#' @return Vector of survival probabilities for the size distribution
calc_surv <- function(z.v, p, n_sz, X.s) {
  z <- z_pow(z.v, n_sz)
  u <- exp(z %*% p$s_z + c(t(X.s) %*% p$s_x))
  u[u<0] <- 0  # in case -Inf
  u[u>1] <- 1  # in case Inf
  return(u / (1+u))
}



##-- Growth
##   z' ~ N(size + environment, sd)
#' Calculate individual growth probabilities based on size & environment
#' @param z1 Vector of sizes to calculate probability of growing TO
#' @param z.v Vector of sizes to calculate probability of growing FROM
#' @param p List of parameters
#' @param n_gz Maximum exponent to raise the size distribution to
#' @param X.g Matrix of environmental covariates 
#' @return Vector of size probabilities for next year for the size distribution
calc_grow <- function(z1, z.v, p, n_gz, X.g) {
  z <- z_pow(z.v, n_gz)
  g <- dnorm(z1, mean=z %*% p$g_z + c(t(X.g) %*% p$g_x), sd=p$g_sig)
  return(g)
}



##-- Flowering
##   P(fl) ~ antilogit(size + environment)
#' Calculate individual flowering probabilities based on size & environment
#' @param z.v Vector of sizes: discretized size range in approximated IPM matrix
#' @param p List of parameters
#' @param n_flz Maximum exponent to raise the size distribution to
#' @param X.fl Matrix of environmental covariates 
#' @return Vector of flowering probabilities for the size distribution
calc_flwr <- function(z.v, p, n_flz, X.fl) {
  z <- z_pow(z.v, n_flz)
  u <- exp(z %*% p$fl_z + c(t(X.fl) %*% p$fl_x))
  u[u<0] <- 0  # in case -Inf
  u[u>1] <- 1  # in case Inf
  return(u / (1+u))
}



##-- Seed production
##   nSeeds ~ exp(size + environment)
#' Calculate individual seed production based on size & environment
#' @param z.v Vector of sizes: discretized size range in approximated IPM matrix
#' @param p List of parameters
#' @param n_seedz Maximum exponent to raise the size distribution to
#' @param X.seed Matrix of environmental covariates 
#' @return Vector of theoretical seed production for the size distribution
calc_seeds <- function(z.v, p, n_seedz, X.seed) {
  z <- z_pow(z.v, n_seedz)
  exp(z %*% p$seed_z + c(t(X.seed) %*% p$seed_x))
}



##-- Recruits (direct)
##   z' ~ P(fl) * nSeeds * P(estab) * N(mn, sd)
#' Calculate size distribution of direct recruits based on size & environment
#' @param z1 Vector of recruit sizes to calculate probability for
#' @param z.v Vector of current population sizes
#' @param p List of parameters
#' @param n_seedz Maximum exponent to raise the size distribution to
#' @param n_flz Maximum exponent to raise the size distribution to
#' @param X.seed Matrix of environmental covariates 
#' @param X.fl Matrix of environmental covariates 
#' @return Vector of recruit size distribution
calc_rcrDir <- function(z1, z.v, p, n_seedz, n_flz, X.seed, X.fl) {
  calc_flwr(z.v, p, n_flz, X.fl) *
    calc_seeds(z.v, p, n_seedz, X.seed) *
    p$rcr_dir *
    dnorm(z1, p$rcr_z[1], p$rcr_z[2])
}



##-- Recruits (seedbank)
##   z' ~ P(rcrSB) * N(mn, sd)
#' Calculate size distribution of seed bank recruits based on size & environment
#' @param z1 Vector of recruit sizes to calculate probability for
#' @param p List of parameters
#' @return Vector of recruit size distribution
calc_rcrSB <- function(z1, p) {
  p$rcr_SB * dnorm(z1, p$rcr_z[1], p$rcr_z[2])
}



##-- Add to seedbank
##   B(z) ~ P(fl) * nSeeds * (1-P(rcrDir)) 
#' Calculate addition to seed bank based on size & environment
#' @param z.v Vector of current population sizes
#' @param p List of parameters
#' @param n_seedz Maximum exponent to raise the size distribution to
#' @param n_flz Maximum exponent to raise the size distribution to
#' @param X.seed Matrix of environmental covariates 
#' @param X.fl Matrix of environmental covariates 
#' @return Vector of seed bank additions by size distribution
calc_addSB <- function(z.v, p, n_seedz, n_flz, X.seed, X.fl) {
  z <- z_pow(z.v, n_seedz)
  calc_flwr(z.v, p, n_flz, X.fl) *
    calc_seeds(z.v, p, n_seedz, X.seed) *
    (1 - p$rcr_dir)
}



##-- Stay in seedbank
##   B ~ P(s.SB) * (1-P(rcrSB))
#' Calculate number of seeds that survive in the seed bank
#' @param p List of parameters
#' @return Probability of failing to recruit & surviving another year in seed bank
calc_staySB <- function(p) {
  p$s_SB * (1 - p$rcr_SB)
}





################################################################################
## IPM matrix construction
########

##-- Set boundaries, meshpoints, & step size for IPM matrix
#' Construct approximate IPM matrix
#' @param n Number of rows and columns in the IPM matrix
#' @param z.rng Vector with two elements defining the typical size range
#' @param buffer Proportion beyond z.rng to allow in individual growth
#' @return List with lo = minimum allowable size, hi = maximum allowable size, 
#' b = matrix boundary points, y = matrix mesh points, h = matrix step size
setup_IPM_matrix <- function(n=100, z.rng=c(1,10), buffer=0) {
  lo <- z.rng[1]*(1-buffer)  # IPM matrix lower limit
  hi <- z.rng[2]*(1+buffer)  # IPM matrix upper limit
  b <- lo + c(0:n)*(hi - lo)/n  # boundary points
  y <- 0.5*(b[1:n] + b[2:(n+1)])  # mesh points
  h <- y[2] - y[1]  # step size
  return(list(lo=lo, hi=hi, b=b, y=y, h=h))
}



##-- Fill P matrix: local 
#' Fill local P matrix, defining size-dependent survival/growth for one grid cell
#' @param h Discretized IPM matrix step size
#' @param y Discretized IPM matrix mesh points
#' @param z.i Vector of sizes: discretized size range in approximated IPM matrix
#' @param p List of parameters
#' @param n_z List of number of size covariates for each vital rate regression
#' @param n_x List of number of environmental covariates for each vital rate 
#' regression
#' @param X_s Matrix of environmental covariates for survival regression
#' @param X_g Matrix of environmental covariates for growth regression
#' @return P matrix with survival & growth kernel
fill_P <- function(h, y, z.i, p, n_z, n_x, X_s, X_g) {
  P.mx <- matrix(0, nrow=p$n+1, ncol=p$n+1)
  # survival & growth
  S.v <- calc_surv(y, p=p, n_sz=n_z$s, X.s=X_s)
  G.mx <- h*outer(y, y, calc_grow, p=p, n_gz=n_z$g, X.g=X_g)
  # correct ejections
  for(k in 1:(p$n/2)) G.mx[1,k] <- G.mx[1,k] + 1 - sum(G.mx[,k])
  for(k in (p$n/2+1):p$n) G.mx[p$n,k] <- G.mx[p$n,k] + 1 - sum(G.mx[,k])
  # fill P matrix
  for(k in z.i) P.mx[k,z.i] <- G.mx[k-1,]*S.v
  P.mx[1,1] <- calc_staySB(p)
  return(P.mx)
}



##-- Fill F matrix: local
#' Fill local F matrix, defining size-dependent fecundity for one grid cell
#' @param h Discretized IPM matrix step size
#' @param y Discretized IPM matrix mesh points
#' @param z.i Vector of sizes: discretized size range in approximated IPM matrix
#' @param p List of parameters
#' @param n_z List of number of size covariates for each vital rate regression
#' @param n_x List of number of environmental covariates for each vital rate regression
#' @param X_fl Matrix of environmental covariates for flowering regression
#' @param X_seed Matrix of environmental covariates for seed regression
#' @param X_germ Matrix of environmental covariates for germination regression
#' @return F matrix with fecundity kernel
fill_F <- function(h, y, z.i, p, n_z, n_x, X_fl, X_seed, X_germ=NULL) {
  F.mx <- matrix(0, nrow=p$n+1, ncol=p$n+1)
  if(!is.null(X_germ)) p$rcr_dir <- p$rcr_SB <- antilogit(c(X_germ %*% p$germ_x))
  F.mx[z.i,z.i] <- outer(y, y, calc_rcrDir, p=p, n_seedz=n_z$seed, n_flz=n_z$fl, 
                         X.seed=X_seed, X.fl=X_fl) * h * p$p_est
  F.mx[1,z.i] <- calc_addSB(y, p=p, n_seedz=n_z$seed, n_flz=n_z$fl, 
                            X.seed=X_seed, X.fl=X_fl) * h
  F.mx[z.i,1] <- calc_rcrSB(y, p) * p$p_est
  return(F.mx)
}



##-- Fill all IPM objects
#' Generates an IPM matrix for each cell, incorporating short distance dispersal
#' if p$p_emig > 0.
#' @param n.cell Number of cells in landscape
#' @param buffer Proportion beyond z.rng to allow in individual growth
#' @param discrete Number of discrete life stages (e.g., seed bank) to include
#' @param p List of parameters
#' @param n_z List of number of size covariates for each vital rate regression
#' @param n_x List of number of environmental covariates for each vital rate regression
#' @param X List of covariates, with elements \code{.$s, .$g, .$fl, .$seed}
#' @param sdd.ji Short distance dispersal immigrant neighborhoods TO each cell j
#' @param p.ji Dispersal probabilities TO each target cell j
#' @param verbose \code{FALSE} Give status updates?
#' @return List of IPMs = IPM matrix for each cell, Ps = P matrix for each cell,
#' Fs = F matrix for each cell, lo = minimum allowable size, hi = maximum 
#' allowable size, b = matrix boundary points, y = matrix mesh points, h = 
#' matrix step size, sdd.j = Short distance dispersal immigrant neighborhoods to 
#' each cell (perspective is the dispersal TO each target cell j)
fill_IPM_matrices <- function(n.cell, buffer, discrete, p, n_z, n_x, 
                              X, sdd.ji, p.ji, verbose=F) {
  library(tidyverse)
  i <- 1:n.cell
  
  # storage objects
  IPMs <- Ps <- Fs <- Fb <- array(0, dim=c(p$n+discrete, p$n+discrete, n.cell))
  
  # IPM matrix setup
  Mx <- setup_IPM_matrix(p$n, p$z.rng, buffer)
  z.i <- (1:p$n) + discrete  # continuous stage matrix indices

  ## local growth
  if(verbose) cat("Beginning local growth \n")
  Ps <- vapply(i, function(x) fill_P(Mx$h, Mx$y, z.i, p, n_z, n_x, 
                                     X$s[x,], X$g[x,]), Ps[,,1])
  Fb <- vapply(i, function(x) fill_F(Mx$h, Mx$y, z.i, p, n_z, n_x, X$fl[x,], 
                                     X$seed[x,], X$germ[x,]), Fb[,,1])
  Fs[z.i,1,] <- Fb[z.i,1,]  # recruits from seedbank unaffected by immigration
  Fs[1,z.i,] <- (1-p$p_emig) * Fb[1,z.i,]  # local contribution to seedbank
  Fs[z.i,z.i,] <- (1-p$p_emig) * Fb[z.i,z.i,]  # local direct recruits
  
  ## dispersal
  if(verbose) cat("Beginning dispersal \n")
  if(p$p_emig > 0) {
    # seed bank: F[1,z,i]
    #   local seeds = Fs[1,z,i]
    #   immigrant seeds = Fb[1,z,j_to_i] * p_emig * pr(j_to_i)
    Fs[1,z.i,] <- vapply(i, function(x) Fs[1,z.i,x] + 
                           Fb[1,z.i,sdd.ji[[x]]] %*% 
                           as.matrix(p$p_emig*p.ji[[x]]),
                         Fs[1,z.i,1])
    # direct recruits: F[z,z,i]
    #   local seedlings = Fs[z,z,i]
    #   immigrant seeds -> seedlings = Fb[z,z,j_to_i] * p_emig * pr(j_to_i)
    D <- lapply(i, function(x) Reduce(`+`, map2(sdd.ji[[x]], p$p_emig*p.ji[[x]], 
                                                ~(Fb[z.i,z.i,.x] * .y))))
    D_NULL <- which(map_lgl(D, is.null))
    if(length(D_NULL) > 0) {
      for(x in D_NULL) D[[x]] <- matrix(0, ncol=p$n, nrow=p$n)
    }
    Fs[z.i,z.i,] <- vapply(i, function (x) Fs[z.i,z.i,x] + D[[x]], Fs[z.i,z.i,1])
  }
  
  return(list(IPMs=Ps+Fs,# Ps=Ps, Fb=Fb, Fs=Fs, 
              lo=Mx$lo, hi=Mx$hi, b=Mx$b, y=Mx$y, h=Mx$h))
}



