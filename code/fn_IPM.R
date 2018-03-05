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
  g <- dnorm(z1, mean=z %*% p$g_z + c(t(X.g) %*% p$g_x), 
             sd=p$g_sig)
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
calc_seeds <- function(z.v, p, n_seeds, X.seed) {
  z <- z_pow(z.v, n_seedz)
  exp(z %*% p$seed_z + c(t(X.seed) %*% p$seed_x))
}


##-- Recruits (direct)
##   z' ~ P(fl) * nSeeds * P(rcrDir) * P(estab) * N(mn, sd)
calc_rcrDir <- function(z1, z.v, p, n_seedz, n_flz, X.seed, X.fl) {
  calc_flwr(z.v, p, n_flz, X.fl) *
    calc_seeds(z.v, p, n_seedz, X.seed) *
    p$rcr_dir *
    p$p_est *
    dnorm(z1, p$rcr_z[1], p$rcr_z[2])
}


##-- Recruits (seedbank)
##   z' ~ P(s.SB) * P(estab) * N(mn, sd)
calc_rcrSB <- function(z1, p) {
  p$s_SB * p$p_est * dnorm(z1, p$rcr_z[1], p$rcr_z[2])
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


