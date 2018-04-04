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
calc_surv <- function(z.v, p, n_sz, X.s) {
  z <- z_pow(z.v, n_sz)
  u <- exp(z %*% p$s_z + c(t(X.s) %*% p$s_x))
  return(u / (1+u))
}



##-- Growth
##   z' ~ N(size + environment, sd)
calc_grow <- function(z1, z.v, p, n_gz, X.g) {
  z <- z_pow(z.v, n_gz)
  g <- dnorm(z1, mean=z %*% p$g_z + c(t(X.g) %*% p$g_x), sd=p$g_sig)
  return(g)
}



##-- Flowering
##   P(fl) ~ antilogit(size + environment)
calc_flwr <- function(z.v, p, n_flz, X.fl) {
  z <- z_pow(z.v, n_flz)
  u <- exp(z %*% p$fl_z + c(t(X.fl) %*% p$fl_x))
  return(u / (1+u))
}



##-- Seed production
##   nSeeds ~ exp(size + environment)
calc_seeds <- function(z.v, p, n_seedz, X.seed) {
  z <- z_pow(z.v, n_seedz)
  exp(z %*% p$seed_z + c(t(X.seed) %*% p$seed_x))
}



##-- Recruits (direct)
##   z' ~ P(fl) * nSeeds * P(estab) * N(mn, sd)
calc_rcrDir <- function(z1, z.v, p, n_seedz, n_flz, X.seed, X.fl) {
  calc_flwr(z.v, p, n_flz, X.fl) *
    calc_seeds(z.v, p, n_seedz, X.seed) *
    p$rcr_dir *
    dnorm(z1, p$rcr_z[1], p$rcr_z[2])
}



##-- Recruits (seedbank)
##   z' ~ P(s.SB) * N(mn, sd)
calc_rcrSB <- function(z1, p) {
  p$s_SB * dnorm(z1, p$rcr_z[1], p$rcr_z[2])
}



##-- Add to seedbank
##   B(z) ~ P(fl) * nSeeds * (1-P(rcrDir)) * P(s.SB)
calc_addSB <- function(z.v, p, n_seedz, n_flz, X.seed, X.fl) {
  z <- z_pow(z.v, n_seedz)
  calc_flwr(z.v, p, n_flz, X.fl) *
    calc_seeds(z.v, p, n_seedz, X.seed) *
    (1 - p$rcr_dir) *
    p$s_SB
}



##-- Stay in seedbank
##   B ~ P(s.SB) * (1-P(rcrDir))
calc_staySB <- function(p) {
  p$s_SB * (1 - p$rcr_SB)
}





################################################################################
## IPM matrix construction
########

##-- Set boundaries, meshpoints, & step size for IPM matrix
setup_IPM_matrix <- function(n=100, z.rng=c(1,10), buffer=0.25) {
  lo <- z.rng[1]*(1-buffer)  # IPM matrix lower limit
  hi <- z.rng[2]*(1+buffer)  # IPM matrix upper limit
  b <- lo + c(0:n)*(hi - lo)/n  # boundary points
  y <- 0.5*(b[1:n] + b[2:(n+1)])  # mesh points
  h <- y[2] - y[1]  # step size
  return(list(lo=lo, hi=hi, b=b, y=y, h=h))
}



##-- Fill P matrix: local 
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



##-- Fill F matrix 
fill_F <- function(h, y, z.i, p, n_z, n_x, X_fl, X_seed) {
  F.mx <- matrix(0, nrow=p$n+1, ncol=p$n+1)
  F.mx[z.i,z.i] <- outer(y, y, calc_rcrDir, p=p, n_seedz=n_z$seed, n_flz=n_z$fl, 
                         X.seed=X_seed, X.fl=X_fl) * h * p$p_est
  F.mx[1,z.i] <- calc_addSB(y, p=p, n_seedz=n_z$seed, n_flz=n_z$fl, 
                            X.seed=X_seed, X.fl=X_fl) * h
  F.mx[z.i,1] <- calc_rcrSB(y, p) * p$p_est
  return(F.mx)
}



##-- Fill all IPM objects
fill_IPM_matrices <- function(n.cell, buffer, discrete, p, n_z, n_x, 
                              X, sdd, sdd.i, lam.final=T, verbose=FALSE) {
  library(tidyverse)
  i <- 1:n.cell
  
  # storage objects
  IPMs <- Ps <- Fs <- Fb <- array(0, dim=c(p$n+1, p$n+1, n.cell))
  Nt <- array(dim=c(p$n+1, n.cell, p$tmax+1))
  sdd.j <- vector("list", length=n.cell)
  lam.t <- p_est.t <- matrix(nrow=n.cell, ncol=p$tmax)
  
  # IPM matrix setup
  L <- setup_IPM_matrix(p$n, p$z.rng, buffer)
  z.i <- (1:p$n) + discrete  # continuous stage matrix indices

  ## local growth
  Ps <- vapply(i, function(x) fill_P(L$h, L$y, z.i, p, n_z, n_x, 
                                            X$s[x,], X$g[x,]), Ps[,,1])
  Fb <- vapply(i, function(x) fill_F(L$h, L$y, z.i, p, n_z, n_x, 
                                            X$fl[x,], X$seed[x,]), Fb[,,1])
  Fs[z.i,z.i,] <- (1-p$p_emig) * Fb[z.i,z.i,]
  Fs[1,z.i,] <- (1-p$p_emig) * Fb[1,z.i,]
  Fs[z.i,1,] <- Fb[z.i,1,]
  if(verbose) cat("Finished local growth \n")
  
  ## dispersal & density dependence
  sdd.j <- lapply(i, function(x) which(sdd[,,2,]==sdd.i[x], arr.ind=T)) 
  p.ij <- lapply(i, function(x) sdd[,,1,][sdd.j[[x]]]) 
  Fs[1,z.i,] <- vapply(i, function(x) Fs[1,z.i,x] + 
                         Fb[1,z.i,sdd.j[[x]][,3]] %*% (p$p_emig*p.ij[[x]]),
                       Fs[1,z.i,1])
  Fs[z.i,z.i,] <- vapply(i, function (x) Fs[z.i,z.i,x] + 
                           Reduce(`+`, map2(sdd.j[[x]][,3], p$p_emig*p.ij[[x]], 
                                            ~(Fb[z.i,z.i,.x] * .y))),
                         Fs[z.i,z.i,1])
  if(verbose) cat("Finished dispersal \n")
  if(p$NDD) {
    Nt[,,1] <- rpois((p$n+1)*n.cell, p$n0/(p$n+1))
    Ft <- Fs/p$p_est  # Fs was already multiplied by p$p_est
    for(k in 1:p$tmax) {
      p_est.t[,k] <- pmin(vapply(i, function(x) p$NDD_n/
                                   sum(Ft[z.i,,x]*Nt[z.i,x,k]), 1), p$p_est)
      Ft[z.i,,] <- vapply(i, function(x) Ft[z.i,,x] * p_est.t[x,k], Ft[z.i,,1])
      IPM.k <- Ps + Ft
      Nt[,,k+1] <- round(sapply(i, function(x) IPM.k[,,x] %*% Nt[,x,k]))
      if(lam.final && k==p$tmax) {
        lam.t[,k] <- vapply(i, function(x) Re(eigen(IPM.k[,,x])$values[1]), 1)
      }
      if(verbose) cat("Finished NDD for year", k, "\n")
    }
  }
  
  return(list(IPMs=Ps+Fs, Ps=Ps, Fb=Fb, Fs=Fs, 
              lo=L$lo, hi=L$hi, b=L$b, y=L$y, h=L$h, sdd.j=sdd.j, 
              Nt=Nt, lam.t=lam.t, p_est.t=p_est.t))
}


