# R code illustrating implementation of an age x size IPM using component
# matrices, as described in Appendix A of Ellner and Rees, Am Nat (2005 or
# thereabouts). The code can be copied and pasted directly into R, or sourced.
#
# The model here is based on our IPM for Onopordum illyricum, but simplified by
# eliminating the variation in individual quality (survival intercept). The
# model is density independent with seedling establishment probability set to
# give a slowly growing population.
#
# Getting on good terms with this model will be a good warmup for working
# through the Onopordum model code. However some of the case-specific tricks
# used to speed up the Onopordum model are not used here, in order to illustrate
# the methods described in Appendix A.

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

# Function to do an image plot of a matrix in the usual orientation, A(1,1) at
# top left
matrix.image <- function(x, y, A, col=topo.colors(100), ...) {
	nx <- length(x); ny <- length(y); 
	x1 <- c(1.5*x[1] - 0.5*x[2], 1.5*x[nx] - 0.5*x[nx-1]); 
	y1 <- c(1.5*y[1] - 0.5*y[2], 1.5*y[ny] - 0.5*y[ny-1]); 
	image(list(x=x, y=y, z=t(A)), xlim=x1, ylim=rev(y1), col=col, bty="u", ...);  
}

#==============================================================================#
# (I) Parameters and demographic functions for computing the kernel 
#==============================================================================#

#probability of seedling establishment
p.est <- 0.03         

# rounded parameter vector, see Table 1 for more details
p.vec <- rep(0,15);
p.vec[1] <- -1.42						#Intercept survival
p.vec[2] <- 1.08						#Size slope survival
p.vec[3] <- -1.09						#Age slope survival
p.vec[4] <- 0.82							#sigma survival: NOT USED HERE  
p.vec[5] <- -24.0						#Intercept flower
p.vec[6] <- 2.9							#Size slope flower
p.vec[7] <- 0.84							#Age slope flower 
p.vec[8] <- 3.2							#Intercept growth
p.vec[9] <- 0.56							#Slope growth
p.vec[10] <- 42.5						#sigma2 growth 
p.vec[11] <- -11.8		  				#Intercept seeds
p.vec[12] <- 2.27		  				#Slope seeds
p.vec[13] <- 1.06						#mean kids size 
p.vec[14] <- 3.37 						#sigma2 kids size 
p.vec[15] <- -0.35					        #varExp para for growth

#Store original p.vec values 
store.p.vec <- p.vec

###### Probability of survival: logistic regression on x=log(size) and age
sx <- function(x, age, params) {
	u <- exp(params[1] + params[2]*x + params[3]*age);
return(u/(1+u));
}

###### Probability of flowering, logistic regression on size and age
fx <- function(x, age, params) {
	u <- exp(params[5] + params[6]*x + params[7]*age);
return(u/(1+u));
}

###### Growth function, normal distribution with mean (mux) and variance
#(sigmax2) depending on size NOTE that in the paper this is called g(y,x). Here
#and in the paper x=size now, y= size next year.
gxy <- function(x, y, age, params) {
	mux <- params[8]+params[9]*x;
	sigmax2 <- params[10]*exp(2*(params[15]*mux));
	sigmax <- sqrt(sigmax2);
	fac1 <- sqrt(2*pi)*sigmax;
	fac2 <- ((y-mux)^2)/(2*sigmax2);
return(exp(-fac2)/fac1);
}

###### Combined survival/growth function, equation 2 in paper 
# NOTE that in the paper this is called p(y,x) 
pxy <- function(x, y, age, params) { 
  return(sx(x, age, params)*(1-fx(x, age, params))*gxy(x, y, age, params)) 
}

# Fecundity function, expected number of recruits of size y produced by a size x
# individual NOTE that in the paper this is called f(y,x)
fxy <- function(x, y, age, params) {
	#expected number of seedlings after establishment
	nkids <- p.est*exp(params[11] + params[12]*x);
	kidsize.mean <- params[13];
	kidsize.var <- params[14]; 
	#probability of producing a seedling of size y
	tmp <- dnorm(y, kidsize.mean, sqrt(kidsize.var))/
	  (1-pnorm(0, kidsize.mean, sqrt(kidsize.var)))
	f <- sx(x, age, params)*fx(x, age, params)*nkids*tmp;
return(f);
}

#==============================================================================#
# (II) Construct the component matrices and their transposes 
#==============================================================================#

# Global variables for midpoint rule approximation 
# n.big.matrix is the number of mesh points for size, n.age is the number of age
# classes.
n.big.matrix <- 100; n.age <- 8;

##### Put all component matrices into 3-dimensional arrays 
#P[j,i,a] will be h*P_{a-1}(x_j,x_i)
P <- array(NA, dim=c(n.big.matrix, n.big.matrix, n.age))
#B[j,i,a] will be h*F_{a-1}(x_j,x_i)
B <- array(NA, dim=c(n.big.matrix, n.big.matrix, n.age)) 

# minimum (0.9*minimum size from data) and maximum sizes (1.1*maximum size from
# data)
minsize <- 0; maxsize <- 9.24; L <- minsize; U <- maxsize;
n <- n.big.matrix 
	
# boundary points b and mesh points y
b <- L+c(0:n)*(U-L)/n; 
y <- 0.5*(b[1:n]+b[2:(n+1)]);

# step size for midpoint rule, see equations 4 and 5
h <- y[2]-y[1]

#Create P and B matrices, looping over age  
for(age in 1:n.age) {
	# Note the transpose is used so we get P(y,x) and B(y,x) from functions that
	# compute p(x,y) and f(x,y). This could have been done better, but we
	# did it this way with Onopordum so you're getting it here. 
	# Note the use of outer to compute each component matrix without looping.    
	P[,,age] <- h * t(outer(y, y, pxy, age=age, params=p.vec))
	B[,,age] <- h * t(outer(y, y, fxy, age=age, params=p.vec))
}

# Create all the transposes
tP <- P; tB <- B; 
for(age in 1:n.age) {
  tP[,,age] <- t(P[,,age]); 
  tB[,,age] <- t(B[,,age]);
} 

#==============================================================================#
# (III) Model iteration functions  
#==============================================================================#
Nt <- matrix(0, n.big.matrix, n.age); 
Nt1 <- Nt; # population now and next year 

iteration <- function(Nt) {
	for(age in 2:n.age) {
	  Nt1[,age] <- P[,,age-1] %*% Nt[,age-1]
	}
	Nt1[,1] <- 0; 
	for(age in 1:n.age) {
	  Nt1[,1] <- Nt1[,1] + B[,,age] %*% Nt[,age]
	}
	return(Nt1)
}

iteration.t <- function(Nt) {
	Nt1[,n.age] <- tB[,,age] %*% Nt[,1]
	for(age in 1:(n.age-1)) {
	  Nt1[,age] <- tB[,,age] %*% Nt[,1] + tP[,,age] %*% Nt[,age+1]
	}
	return(Nt1)
} 

#==============================================================================#
# (IV) Start using the model 
#==============================================================================#

##### Estimate lambda and w by iterating unperturbed matrix  
Nt <- matrix(1, n.big.matrix, n.age); 
qmax <- 1000; lam <- 1; tol <- 1.e-8; 
while(qmax>tol) {
	Nt1 <- iteration(Nt);
	qmax <- sum(abs(Nt1-lam*Nt));  
	lam <- sum(Nt1); 
	Nt <- Nt1/lam; 
	cat(lam, qmax, "\n");
} 
stable.dist <- Nt/sum(Nt); 
lam.stable <- lam; 	

##### Estimate lambda and v by transpose iteration 
Nt <- matrix(1, n.big.matrix, n.age); 
qmax <- 1000; lam <- 1; tol <- 1.e-8; 
while(qmax>tol) {
	Nt1 <- iteration.t(Nt);
	qmax <- sum(abs(Nt1-lam*Nt));  
	lam <- sum(Nt1); 
	Nt <- Nt1/lam; 
	cat(lam, qmax, "\n");
} 
repro.val <- Nt/sum(Nt); 
lam.stable.t <- lam; 

#Check that forward and transpose iteration give the same lambda 
cat("Forward lambda:", lam.stable,
    "should equal transpose lambda:",lam.stable.t,"\n"); 	

####### Sensitivities/elasticities from sensitivity formula 
## Note: s(y,x) is computed only for transitions with K(y,x)>0, hence e(y,x)>0.  

# Survival-growth transitions go to age=(current age+1), all sizes   
Psens=array(0,dim=c(n.big.matrix,n.big.matrix,n.age))
# Psens[,,j] will hold the s(y,x) values for size transitions between ages j-1 and j. 
for(age in 1:(n.age-1)) {
	Psens[,,age]=outer(repro.val[,age+1],stable.dist[,age])
}
v.dot.w = h*sum(repro.val*stable.dist); 
Psens=Psens/v.dot.w; 
Pelas=Psens*(P/h)/lam; 
rm(Psens) # it's big and not needed any more, so get it out of memory.   

# Fecundity transitions go from all (size x age) combinations to all sizes at age=0. 
Bsens=array(0,dim=c(n.big.matrix,n.big.matrix,n.age))
# Bsens[,,j] will hold the s(y,x) values for size-specific fecundities of age j-1 parents. 
	for(age in 1:n.age) {
		Bsens[,,age]=outer(repro.val[,1],stable.dist[,age]); 
	}
	
Bsens=Bsens/v.dot.w; 
Belas=Bsens*(B/h)/lam; 
rm(Bsens);  

#check elasticities integrate to 1 as they should (by midpoint rule) 
sum(h*h*Pelas)+sum(h*h*Belas);
cat("Integral of elasticity function = ",sum(h*h*Pelas)+sum(h*h*Belas),"\n"); 

#============================================================================# 
#  Plot stable distribution and size-dependent total elasticity 
#============================================================================# 
par(mfrow=c(2,2)); 

plot(y,apply(stable.dist,1,sum),xlab="Size x",ylab="frequency",type="l");
title(main="Stable size distribution");  
barplot(apply(stable.dist,2,sum),width=1,names.arg=as.character(0:7),
xlab="Age a",ylab="frequency",type="l"); 
title(main="Stable age distribution"); 

Pelas.size=apply(Pelas,c(1,2),sum); #sum elasticities over age  
matrix.image(y,y,Pelas.size,xlab="Size t",ylab="Size t+1");
contour(y,y,t(Pelas.size),add=T,cex=3);
title(main="Survival-growth elasticities: 74%"); 

Belas.size=apply(Belas,c(1,2),sum); #sum elasticities over age 
matrix.image(y,y,Belas.size,xlab="Size t",ylab="Size t+1");
contour(y,y,t(Belas.size),add=T,cex=3);
title(main="Fecundity elasticities: 26%"); 

