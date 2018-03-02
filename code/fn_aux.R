# Auxilliary functions
# Comparison of SDM approaches
# Tim Szewczyk


##-- set boundaries, meshpoints, step size for IPM matrix
##   create objects in global environment
setup_IPM_matrix <- function(n=100, z.rng=c(1,10), buffer=0.25, discrete=0) {
  lo <<- z.rng[1]*(1-buffer)  # IPM matrix lower limit
  hi <<- z.rng[2]*(1+buffer)  # IPM matrix upper limit
  b <<- lo + c(0:n)*(hi - lo)/n  # boundary points
  y <<- 0.5*(b[1:n] + b[2:(n+1)])  # mesh points
  h <<- y[2] - y[1]  # step size
  z.i <<- (1:n) + discrete  # continuous stage matrix indices
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
