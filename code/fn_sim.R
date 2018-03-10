



##--- simulate expected values
## 
sim_expected <- function(d_i, p, X, n_z) {
  E <- list(s=NULL, g=NULL, fl=NULL, seed=NULL)
  size.mx <- cbind(1, d_i$size, d_i$size^2, d_i$size^3)
  E$s <- antilogit(size.mx[,1:n_z$s] %*% p$s_z + c(X$s %*% p$s_x))
  E$g <- size.mx[,1:n_z$g] %*% p$g_z + c(X$g %*% p$g_x)
  E$fl <- antilogit(size.mx[,1:n_z$fl] %*% p$fl_z + c(X$fl %*% p$fl_x))
  E$seed <- exp(size.mx[,1:n_z$seed] %*% p$seed_z + c(X$seed %*% p$seed_x))
  return(E)
}


##--- simulate realized values
##
sim_realized <- function(d_i, k.index, E, p, lo, hi) {
  n_t <- length(k.index)
  d_i$surv[k.index] <- rbinom(n_t, 1, E$s)
  d_i$sizeNext[k.index] <- rnorm(n_t, E$g, p$g_sig) * 
    ifelse(d_i$surv[k.index], 1, NA)
  d_i$sizeNext[k.index] <- ifelse(d_i$sizeNext[k.index] < lo, 
                                  lo, d_i$sizeNext[k.index])
  d_i$sizeNext[k.index] <- ifelse(d_i$sizeNext[k.index] > hi, 
                                  hi, d_i$sizeNext[k.index])
  d_i$fl[k.index] <- rbinom(n_t, 1, E$fl) * ifelse(d_i$surv[k.index], 1, NA)
  d_i$seed[k.index] <- rpois(n_t, E$seed) * 
    ifelse(d_i$surv[k.index] & d_i$fl[k.index], 1, NA)
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

