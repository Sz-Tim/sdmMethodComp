# R code illustrating implementation of a size x quality IPM as described in
# Appendix A of Ellner and Rees, Am Nat (2005 or thereabouts). 
# The code can be copied and pasted directly into R, or sourced.
#
# The model here starts from the size-structured model in CompareSensitivity.R 
# and adds individual quality as a second dynamic i-state variable. The quality
# dynamics are a linear autoregressive process, as in the Pfister-Stevens model
# described in section 3 of Appendix A, and we use the 'stacking' procedure
# suggested there to create a single large iteration matrix. You then have to
# be careful about 'unstacking' to plot results and do sensitivity analysis. 
#
# This code works with version 0.98-7 of the Matrix package; with 
# earlier versions it might fail (later ones, too...) 

#clear everything
rm(list=ls(all=TRUE))

require(Matrix); 

# Set matrix size (to show up errors) and convergence tolerance. 
# Make m1 and m2 smaller to make this run faster.  
m1=51; m2=49; tol=1.e-8; 

# Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
matrix.image=function(x,y,A,col=topo.colors(200),...) {
	nx=length(x); ny=length(y); 
	x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
	y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
	image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}


#============================================================================# 
#  Define the kernel K(y,x) and iteration matrix, for x=size in range [-2,8]
#============================================================================# 
# Survival probability is logistic regression model, s(x)=exp(x-1)/1+exp(x-1)
# New size = Normal(mean=1+0.7*old size+quality,sd=0.3)
# New quality = Normal(mean=.5*old quality,sd=.5)  
# Fecundity = 0.75 per surviving individual
# Offspring size = Normal(mean=0.8+0.2*parent size, sd=0.35)
# Offspring quality=Normal(mean=0,sd=.5)
      
pyx=function(xp,qp,x,q) {exp(x-1)/(1+exp(x-1))*dnorm(xp,mean=1+0.7*x+q,sd=0.3)*dnorm(qp,mean=0.5*q,sd=0.5)}
fyx=function(xp,qp,x,q) {0.75*exp(x-1)/(1+exp(x-1))*dnorm(xp,mean=0.8+x/5,sd=0.35)*dnorm(qp,mean=0,sd=0.5)}
kyx=function(xp,qp,x,q) {pyx(xp,qp,x,q)+fyx(xp,qp,x,q)}

# Compute meshpoints  
h1=9/m1; y1=(h1/2)*((0:(m1-1))+(1:m1))-2;
h2=5/m2; y2=(h2/2)*((0:(m2-1))+(1:m2))-2.5;

# Compute the iteration matrix. With a bit of vectorizing it's not too slow,
# though you can probably do better if you need to. The shortcuts here have 
# been checked against the results from code that uses loops for everything. 
plop=function(i,j) {(j-1)*m1+i} # for putting values in proper place in A 
Plop=outer(1:m1,1:m2,plop); 

A=matrix(0,m1*m2,m1*m2); Kvals=array(0,c(m1,m2,m1,m2));  
for(i in 1:m1){
	for(j in 1:m2){
		for(k in 1:m1){
				kvals=kyx(y1[k],y2[1:m2],y1[i],y2[j])
				A[Plop[k,1:m2],Plop[i,j]]=kvals
				Kvals[k,1:m2,i,j]=kvals
			
	}}
	cat(i,"\n"); 
}
h12=h1*h2; A=h12*A;  

#============================================================================# 
#  Find lambda, w by iteration 
#
#  Note: the Matrix package is used to iterate more quickly. The check  
#  for convergence requires extracting matrix entries via the @x slot
#  of a Matrix object. Matrix is S4-style -- see ?Matrix. 
#============================================================================# 
A2=Matrix(A); nt=Matrix(1,m1*m2,1); nt1=nt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	nt1=A2%*%nt;
	qmax=sum(abs((nt1-lam*nt)@x));  
	lam=sum(nt1@x); 
	nt@x=(nt1@x)/lam; #we're cheating here - don't tell Doug Bates.  
	cat(lam,qmax,"\n");
} 
nt=matrix(nt@x,m1*m2,1); 
stable.dist=nt/(h1*h2*sum(nt)); #normalize so that integral=1
lam.stable=lam; 

# Check that the @bits worked as intended.   
qmax=sum(abs(lam*nt-A%*%nt)); 
cat("Convergence: ",qmax," should be less than ",tol,"\n"); 

#============================================================================# 
#  Find reproductive value function by iteration.  
#  Note that we use left-multiplication to do the transpose iteration,
#  which is equivalent to using the transpose kernel. 
#============================================================================# 
vt=Matrix(1,1,m1*m2); vt1=vt; 

qmax=1000; lam=1; 
while(qmax>tol) {
	vt1=vt%*%A2;
	qmax=sum(abs((vt1-lam*vt)@x));  
	lam=sum(vt1@x); 
	vt@x=(vt1@x)/lam;   
	cat(lam,qmax,"\n");
} 
v=t(matrix(vt@x,1,m1*m2)); 
lam.stable.t=lam; 

# Compute elasticity matrix 
repro.val=matrix(v,m1,m2); stable.state=matrix(stable.dist,m1,m2); 
v.dot.w=sum(h1*h2*stable.state*repro.val)
sens=outer(repro.val,stable.state)/v.dot.w
elas=sens*Kvals/lam.stable;

# Compute matrix of total(=integrated) elasticities for all transitions (x_i,q_j) -> anywhere 
total.elas=h1*h2*apply(elas,c(3,4),sum); 

# Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity=",sum(h1*h2*h1*h2*elas)," Should =1","\n"); 

#============================================================================# 
#  Plot stable state distribution and state-dependent total elasticity 
#============================================================================# 
par(mfrow=c(2,2)); 

plot(y1,apply(stable.state,1,sum),xlab="Size x",ylab="frequency",type="l"); 
title(main="Stable size distribution");

plot(y2,apply(stable.state,2,sum),xlab="Quality q",ylab="frequency",type="l"); 
title(main="Stable quality distribution"); 

matrix.image(y2,y1,stable.state,xlab="Quality q",ylab="Size x");
contour(y2,y1,t(stable.state),add=T,cex=3);
title(main="Stable size-quality distribution"); 

matrix.image(y2,y1,total.elas,xlab="Quality q",ylab="Size x");
contour(y2,y1,t(total.elas),add=T,cex=3);
title(main="State-dependent total elasticity"); 


