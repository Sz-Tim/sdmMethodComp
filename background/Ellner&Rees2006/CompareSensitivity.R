# R code illustrating sensitivity analysis by perturbation 
# from Ellner and Rees, Am Nat (2005 or thereabouts). 
# See Appendix A section 5 of the paper for more details.
# The code can be copied and pasted directly into R, or sourced.
#
# For purposes of illustration we use a basic IPM with continuous size structure. 
# The iteration matrix is relatively small (bigM x bigM) so everything goes fast, 
# and we can use a very strict convergence tolerance (tol=1.e-10) on iterations  
# in order to get very accurate results. In practice with larger iteration
# matrices you could start with a fairly loose tol, and decrease it until results
# stabilize.

#Set matrix size and convergence tolerance 
bigM=250; tol=1.e-10 

# Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
matrix.image=function(x,y,A,col=topo.colors(100),...) {
	nx=length(x); ny=length(y); 
	x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
	y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
	image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}

#============================================================================# 
#  Define the kernel K(y,x) and iteration matrix, for x=size in range [0,5]
#============================================================================# 
# Survival probability is logistic regression model, s(x)=exp(x-1)/(1+exp(x-1))
# New size = Normal(mean=1+0.7*old size,sd=0.3)
# Fecundity = 0.75 per surviving individual
# Offspring size = Normal(mean=0.8+0.2*parent size, sd=0.35)      
pyx=function(y,x) {exp(x-1)/(1+exp(x-1))*dnorm(y,mean=1+0.7*x,sd=0.3)}
fyx=function(y,x) {0.75*exp(x-1)/(1+exp(x-1))*dnorm(y,mean=0.8+x/5,sd=0.35)}
kyx=function(y,x) {pyx(y,x)+fyx(y,x)}

# Compute meshpoints iteration matrix KD 
# Note the use of outer() to compute kernel values at all meshpoints in one statement.  
h=5/bigM; y=(h/2)*((0:(bigM-1))+(1:bigM));  
K=outer(y,y,kyx); KD=h*K;  

# Plot the kernel and overlay contours 
par(mfrow=c(2,2),cex.axis=1.35,cex.lab=1.35,bty="o"); 
matrix.image(y,y,K, xlab="year t",ylab="year t+1",col=topo.colors(100)); title("Kernel"); 
contour(y,y,t(K),add=T,cex=3);
abline(h=0); abline(v=0); 

#============================================================================# 
# Sensitivity analysis I: from the formulas 
#============================================================================# 
# Get lamda,v,w from the iteration matrix, and plot 
lam.eigen=as.real(eigen(KD)$values[1]); 
w.eigen=as.real(eigen(KD)$vectors[,1]); w.eigen=w.eigen/sum(w.eigen); 
v.eigen=as.real(eigen(t(KD))$vectors[,1]); v.eigen=v.eigen/v.eigen[1]; 
plot(y,w.eigen,type="l",xlab="Size x",ylab="w (solid), v (dash)"); title("Eigenvectors w and v");
points(y,v.eigen*mean(w.eigen)/mean(v.eigen),type="l",lty=2);  

# Compute sensitivity and elasticity using sensitivity formulas 
v.dot.w= h*sum(v.eigen*w.eigen) # note <v,w> is an integral, done here by midpoint rule.  
sens.eigen=outer(v.eigen,w.eigen)/v.dot.w;   
elas.eigen=K*sens.eigen/lam.eigen; 

#============================================================================# 
# Sensitivity analysis II: do it all by iteration 
#============================================================================# 
# Estimate lambda and w by iterating unperturbed matrix  
nt=matrix(1,bigM); 
qmax=1000; lam=1; 
while(qmax>tol) {
	nt1=KD%*%nt;
	qmax=sum(abs(nt1-lam*nt));  
	lam=sum(nt1); 
	nt=nt1/lam; 
	cat(lam,qmax,"\n");
} 
stable.dist=nt/sum(nt); lam.stable=lam; 

# Sensitivity and elasticity by perturbing the kernel
# Note: entry (i,j) of iteration matrix KD is h*K(y_i,x_j), so each row of KD corresponds
# to kernel values for a particular value of y  
sen.big<-array(NA,dim=c(bigM,bigM)) 	#array to store the results 
for(row in 1:bigM) {			# loop over y values 
	# choose x* to maximize e(y,x) for this y value, by scanning across the row    
	big.one=which(KD[row,]*stable.dist==max(KD[row,]*stable.dist)); 

	# perturb the kernel up and down near (y,x*)
	delta=0.1*h*KD[row,big.one];
	KDup=KD; KDup[row,big.one] = KD[row,big.one]+delta/h;
	KDdown=KD; KDdown[row,big.one] = KD[row,big.one]-delta/h;  
	##### Why delta over h? 
	##### K is perturbed by a function f=delta/(h*h) on grid cell ij,
	##### whose integral over the grid cell is therefore delta.  
	##### Perturbing K by f, perturbs KD by h*f = delta/h. 

	# Now compute lambda for the up- and down-perturbed kernels 
	# As suggested in the paper, start at the unperturbed stable distribution
	# to get rapid convergence to perturbed stable distribution and growth rate 
	qmax=1; lamup=1; lamdown=1; 
	Nt.up<-stable.dist; Nt.down<-stable.dist
	while (qmax>tol){
		Nt1.up=KDup%*%Nt.up
		qmax=sum(abs(Nt1.up-lamup*Nt.up))
		lamup=sum(Nt1.up)
		Nt.up=Nt1.up/lamup
		
		Nt1.down=KDdown%*%Nt.down
		qmax=qmax+sum(abs(Nt1.down-lamdown*Nt.down))
		lamdown<-sum(Nt1.down)
		Nt.down<-Nt1.down/lamdown
	}
	sen.big.row<-(lamup-lamdown)/(2*delta) #sensitivity for perturbation at (y,x*)
	sen.big[row,]<-(stable.dist/stable.dist[big.one])*sen.big.row #sensitivity at other x's 
	#cat(row,big.one," sens=",sen.big.row, "\n")
}
elas.big=K*sen.big/lam.stable;

#============================================================================# 
# Plot the results 
#============================================================================# 
matrix.image(y,y,elas.eigen,xlab="year t",ylab="year t+1",col=topo.colors(100)); 
title("Elasticity from eigenvectors"); 
contour(y,y,t(elas.eigen),add=T,cex=3);
abline(h=0); abline(v=0); 

matrix.image(y,y,elas.big,xlab="year t",ylab="year t+1",col=topo.colors(100)); 
title("Elasticity from iteration"); 
contour(y,y,t(elas.big),add=T,cex=3);
abline(h=0); abline(v=0); 


#============================================================================# 
# Print some checks to the console window
#============================================================================#  

# Compute maximum relative error in elasticity 
relerr=(elas.big-elas.eigen)/elas.eigen; 

# Check that elasticity function by perturbation integrates to 1
elas.big=K*sen.big/lam.stable;

cat("Integrated elasticity from formulas=", h*h*sum(elas.eigen), ", which should be 1", "\n");  
cat("Integrated elasticity by perturbation=", h*h*sum(elas.big), ", which should be 1", "\n");  
cat("Maximum relative error in elasticities is ",round(100*max(relerr),3), " percent", "\n"); 



