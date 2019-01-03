# R code for implementing the Onoporum model from Ellner and Rees, Am Nat (2005 or thereabouts). 
# See Appendix A section 4 of the paper for more details.
# The code can be copied and pasted directly into R, or sourced.

# The state of the population is stored in two three dimensional arrays, 
#	Nt="this year" and Nt1="next year" 
# The indices on these arrays are size, age and quality, respectively.
# Note that age indices 1 to 8 (called "age" in the code) correspond to actual ages 0 to 7.  

### There are 5  sections to the code, headed by #=========# lines . 
# (I)   Parameter values 
# (II)  Utility functions to construct the kernel 
# (III) Functions to iterate the model and do transpose iteration for reproductive value  
# (IV)  Construction of the matrices used to iterate the model 
# (V)   The calculations for figures 1-4 in the text 

### So that it runs quickly, this code uses very few mesh points for size and quality
### For accurate results we recommend using 
	# n.big.matrix <- 50; n.sigma <- 40; 
### You may need to add --max-mem-size=2G to the end of the target line of your R shortcut 
### to use this number of mesh points.

#============================================================================================#
# (I) Parameter Values 
#============================================================================================#
#clear everything
rm(list=ls(all=TRUE))

# Global variables for midpoint rule approximation 
# n.big.matrix and n.sigma are the number mesh points for size and quality, respectively.
# n.age is the number of age classes. 
n.big.matrix <- 25; n.sigma <- 20; n.age <- 8;

# minimum (0.9*minimum size from data) and maximum sizes (1.1*maximum size from data)
minsize=0; maxsize=9.24;

#probability of seedling establishment assuming lambda is roughly 1.026 as in data
p.est<-0.025         

# Rounded parameter vector, see Table 1 for more details. 
p.vec<-rep(0,15);
p.vec[1]<-  -1.42						#Intercept survival
p.vec[2]<-  1.08						#Size slope survival
p.vec[3]<- -1.09						#Age slope survival
p.vec[4]<- 0.82							#sigma survival 
p.vec[5]<- -24.0						#Intercept flower
p.vec[6]<- 2.9							#Size slope flower
p.vec[7]<- 0.84							#Age slope flower 
p.vec[8]<-  3.2							#Intercept growth
p.vec[9]<- 0.56							#Slope growth
p.vec[10]<- 42.5						#sigma2 growth 
p.vec[11]<- -11.8		  				#Intercept seeds
p.vec[12]<- 2.27		  				#Slope seeds
p.vec[13]<- 1.06						#mean kids size 
p.vec[14]<- 3.37 						#sigma2 kids size 
p.vec[15]<- -0.35					        #varExp para for growth

#Store original p.vec values 
store.p.vec<-p.vec

#============================================================================================#
# (II) Utility functions to construct the kernel  
#============================================================================================#
## Functions to compute kernel components from the fitted models

# Probability of survival: logistic regression on x=log(size) and age
sx<-function(x,age,params) {
	u<-exp(params[1]+params[2]*x+params[3]*age);
return(u/(1+u));
}

# Probability of flowering, logistic regression on size and age
fx<-function(x,age,params) {
	u<-exp(params[5]+params[6]*x+params[7]*age);
return(u/(1+u));
}

# Growth function, normal distribution with mean (mux) and variance (sigmax2) depending on size
# NOTE that in the paper this is called g(y,x). Here and in the paper x=size now, y= size next year. 
gxy<-function(x,y,age,params) {
	mux<-params[8]+params[9]*x;
	sigmax2<-params[10]*exp(2*(params[15]*mux));
	sigmax<-sqrt(sigmax2);
	fac1<-sqrt(2*pi)*sigmax;
	fac2<-((y-mux)^2)/(2*sigmax2);
return(exp(-fac2)/fac1);
}

# Combined survival/growth function, equation 2 in paper 
# NOTE that in the paper this is called p(y,x) 
pxy<-function(x,y,age,params) { return(sx(x,age,params)*(1-fx(x,age,params))*gxy(x,y,age,params)) }

# Fecundity function, expected number of recruits of size y produced by a size x individual
# NOTE that in the paper this is called f(y,x) 
fxy<-function(x,y,age,params) {
	#expected number of seedlings after establishment
	nkids<-p.est*exp(params[11]+params[12]*x);
	kidsize.mean<- params[13];
	kidsize.var<- params[14]; 
	#probability of producing a seedling of size y
	tmp<-dnorm(y,kidsize.mean,sqrt(kidsize.var))/(1-pnorm(0,kidsize.mean,sqrt(kidsize.var)))
	f<-sx(x,age,params)*fx(x,age,params)*nkids*tmp;
return(f);
}

#============================================================================================#
# (III) Model iteration functions  
#============================================================================================#
### Iterate the model one time step using B or B.short (created below), starting
### from current population array Nt. 
### Probability of seedling establishment can be density dependent (DI=F) or density independent (DI=T)
iteration<-function(Nt,use.B.short=T,DI=T){
	Nt1<-array(0,dim=c(n.big.matrix,n.age,n.sigma))

	#calculate all births and put in appropriate places
	if (use.B.short==T){
		if (DI==F){
		#density dependent model. The fecundity function fxy includes the density independent p.est, so divide through by p.est
		#to get the total number of seedling *before* density dependence
		total.seeds<-sum(B.short * Nt)/p.est
		#91 is the average number of recruits into the population
		DD.p.est<<-91/total.seeds
		#check if probability of establishment is greater than 1, if so make it 1
		if(DD.p.est>1) DD.p.est<-1
		total.rec<-DD.p.est*total.seeds
		} else {
		total.rec<-sum(B.short * Nt)
		}
		#put new recruits into appropriate places using dist.sigma.d.by.offspring.size to allocate seedling to recruit size by 
		#quality distribution
		Nt1[,1,]<-dist.sigma.d.by.offspring.size*total.rec
	} else {
		# not using B.short, so calculate the number of recruits for each age x quality combination and then allocate these across
		# the recruit size and quality distribution. 
		# Note: it is much quicker to use B.short, this is really just a check it works
		for (inter.p in 1:n.sigma){
 		for (age.p in 1:n.age){
			biths.age.sigma.p<-B[,,age.p,inter.p] %*% Nt[,age.p,inter.p]
			for (inter.b in 1:n.sigma){
				Nt1[,1,inter.b]<-Nt1[,1,inter.b]+prob.sigma.d[inter.b]*biths.age.sigma.p
			}
		}
		}
	}

	# Death and aging. This uses equation 2 to calculate the size distribution of age i individuals for 
	# each quality class (survival intercept). 
	for (inter in 1:n.sigma){
	for (age in 2:n.age){
		Nt1[,age,inter]<-P[,,age-1,inter] %*% Nt[,age-1,inter]
	}
	}
return(Nt1)
}

# Iterate the transpose using formulae from Appendix A section 3
iteration.t<-function(Nt){
	Nt1<-array(0,dim=c(n.big.matrix,n.age,n.sigma))

	for (inter in 1:n.sigma){
	for (age in 1:n.age){
	for (inter.s in 1:n.sigma){
		Nt1[,age,inter]<-Nt1[,age,inter]+prob.sigma.d[inter.s]*t(B[,,age,inter]) %*% Nt[,1,inter.s]
	}
	}
	}

	for (inter in 1:n.sigma){
	for (age in 1:(n.age-1)){ 
		Nt1[,age,inter]<-Nt1[,age,inter]+t(P[,,age,inter]) %*% Nt[,age+1,inter]
	}
	}

return(Nt1)
}

#============================================================================================#
# (IV) Set up Matrices for model iteration 
# P=survival-growth, B=fecundity, each of these is a 4 dimensional array
# with indices for size next year and size, age, quality this year  
#============================================================================================#

P<-array(NA,dim=c(n.big.matrix,n.big.matrix,n.age,n.sigma))
B<-array(NA,dim=c(n.big.matrix,n.big.matrix,n.age,n.sigma))

L<-minsize; U<-maxsize;
n<-n.big.matrix
	
# boundary points b and mesh points y
b<-L+c(0:n)*(U-L)/n;
y<-0.5*(b[1:n]+b[2:(n+1)]);

# step size for mid point rule, see equations 4 and 5
h<-y[2]-y[1]

#distribution of offspring sizes: normal truncated at minsize=0. 
prob.offspring<-dnorm(y,p.vec[13],sqrt(p.vec[14]))
#normalised to sum to 1 
prob.offspring<-prob.offspring/sum(prob.offspring)

#points for numerical integration of survival intercepts
cuts.sigma.d<-seq(-4*p.vec[4],4*p.vec[4],length=n.sigma);
#normalised to sum to 1
prob.sigma.d<-dnorm(cuts.sigma.d)/sum(dnorm(cuts.sigma.d));

#Create P and B matrices, looping over quality class (survival intercept) and age  
for (inter in 1:n.sigma){
 for (age in 1:n.age){
	# Calculate survival intercept for this quality class
	p.vec[1]<-store.p.vec[1]+cuts.sigma.d[inter]
	# Calculate survial and birth matrices for each quality (intercept) age combination. 
	# Note the transpose is used so we get P(y,x) and B(y,x) from functions that
	# compute p(x,y) and f(x,y). This could have been done better, but we're
	# giving you the code that we actually used. 
	P[,,age,inter]<-(U-L)*t(outer(y,y,pxy,age=age,params=p.vec))/n
	B[,,age,inter]<-(U-L)*t(outer(y,y,fxy,age=age,params=p.vec))/n
}
}

# Create B.short - the total number of offspring produced by an individual of a given size, age and quality. 
# The sum of the element wise product (B.short*Nt) gives total number of seeds produced, see Appendix A section 4
B.short<-apply(B,c(2,3,4),sum)

# matrix with elements describing the (size x quality) probability distribution of offspring. 
# Size and quality are assumed to be independent
dist.sigma.d.by.offspring.size<-outer(prob.offspring,prob.sigma.d)	

#============================================================================================#
# (V) Working with the model  
#============================================================================================#

################ iterate model using B.short to get lam and stable.dist 
Nt<-array(1,dim=c(n.big.matrix,n.age,n.sigma));
Nt=Nt/sum(Nt); 

qmax=1; lam=1; iter=1;
tol=1e-6; # convergence tolerance for iteration to stable distribution  

while (qmax>tol){
	Nt1=iteration(Nt,use.B.short=T,DI=T)
	#calculate sum of absolute deviations between Nt1 and lam*Nt - used to measure convergence to stable state
	qmax=sum(abs(Nt1-lam*Nt)); 
	lam=sum(Nt1); 
	cat(iter,sum(Nt),"  ",sum(Nt1),"  ",lam, qmax,"\n")
	iter=iter+1; 
	#set Nt equal to NT1 and normalise so Nt sums to 1
	Nt = Nt1/lam;
}
lam.stable=lam; stable.dist=Nt; 

win.graph()
par(mfrow=c(2,2),bty="l",pch=19)
plot(y,apply(Nt,1,sum),xlab="Size",ylab="Stable size distribution",type="l")
plot(1:n.age,apply(Nt,2,sum),xlab="Age",ylab="Stable age distribution",type="l")
plot(store.p.vec[1]+cuts.sigma.d,apply(Nt,3,sum),xlab="Individual intercept",ylab="Stable intercept distribution",type="l")

#add initial intercept distribution to the plot 
points(store.p.vec[1]+cuts.sigma.d,prob.sigma.d,type="l",col="red")

################ Iterate transpose to get state-dependent reproductive value function v(x) 
# There's no 'big B' trick to use here, so this is slow. 
vt<-array(1,dim=c(n.big.matrix,n.age,n.sigma));
vt=vt/sum(vt); 

qmax=1; lam=1; iter=1;
tol=1e-6; 

while (qmax>tol){
	vt1=iteration.t(vt)
	qmax=sum(abs(vt1-lam*vt)); 
	lam=sum(vt1);
	cat(iter,sum(vt),"  ",sum(vt1),"  ",lam, qmax,"\n")
	iter=iter+1; 
	vt = vt1/lam;
}
lam.t=lam; repro.val=vt;  

rm(vt,vt1)

# Plot stable size distribution and size dependent reproductive value
win.graph()
par(mfrow=c(1,2),bty="l",pty="s")
N.size<-apply(Nt,1,sum)

plot(y,N.size,xlab="Rosette area (log scale)",type="l",ylab="Proportion")
plot(y,apply(repro.val,1,sum), type="l",xlab="Rosette area (log scale)",ylab="Reproductive value",log="y")


#================================================================================
#####    Sensitivity and Elasticity calculations using the sensitivity formula 
#####
##### Note s(y,x) is computed only for transitions with K(y,x)>0, hence e(y,x)>0. 
##### The use of outer() to compute many s(y,x) at once is much faster than looping. 
#####
##### Because of age structure each possible transition is either due to
##### survival-growth only or due to fecundity only, rather than some of each. 
##### We can therefore compute s(y,x) for each part separately using the 
##### sensitivity formula (13): fecundity if age=0 in y, survival/growth if age>0.
#=================================================================================

#################### Survival-growth sensitivities/elasticities
### Survival-growth transitions go to current age+1, current quality class, all sizes   
Psens=array(0,dim=c(n.big.matrix,n.big.matrix,n.age,n.sigma))
for(age in 1:(n.age-1)) {
	for(inter in 1:n.sigma) {
		Psens[,,age,inter]=outer(repro.val[,age+1,inter],stable.dist[,age,inter])
}
}
v.dot.w = h*sum(repro.val*stable.dist); 
Psens=Psens/v.dot.w; 
Pelas=Psens*(P/h)/lam; 
rm(Psens) # it's big, get it out of memory now that it's not needed. 

################### Fecundity sensitivities/elasticities
### Fecundity transitions go to all qualities and all sizes at age=0. 
BigBsens=array(0,dim=c(n.sigma,n.big.matrix,n.big.matrix,n.age,n.sigma))
	for(age in 1:n.age) { #loop on age of parent
	for(inter.p in 1:n.sigma) {#loop on quality of parent
	for(inter.b in 1:n.sigma) {#loop on quality of offspring 
		BigBsens[inter.b,,,age,inter.p]=outer(repro.val[,1,inter.b],stable.dist[,age,inter.p]); 
	}
	}
	}
BigBsens=BigBsens/v.dot.w; 

BigBelas=array(0,dim=c(n.sigma,n.big.matrix,n.big.matrix,n.age,n.sigma))
for(inter.b in 1:n.sigma) {
	BigBelas[inter.b,,,,]=BigBsens[inter.b,,,,]*(prob.sigma.d[inter.b]*B/h)/lam;
}

rm(BigBsens); gc(); # clear more big stuff

#check elasticities integrate to 1 as they should 
sum(h*h*Pelas)+sum(h*h*BigBelas);

# plots for P elasticities
# Note use of locator() to add text (commented out for convenience). 
win.graph()
par(mfrow=c(2,2),pty="s",bty="l",cex.axis=1.35,cex.lab=1.35)
elas.sum<-apply(Pelas,3,sum)*h*h	
plot(0:7,elas.sum,type="h",lwd=10,xlab="Age (years)",ylab="Elasticity")
#text(locator(1),"a)",cex=1.35)
elas.sum<-apply(Pelas,4,sum)*h*h
plot(store.p.vec[1]+cuts.sigma.d,elas.sum,type="h",lwd=2,xlab="Survival intercept",ylab="Elasticity")
#text(locator(1),"b)",cex=1.35)

# can add intercept distribution in offspring to the plot 
# points(store.p.vec[1]+cuts.sigma.d,prob.sigma.d,type="l",col="red")
px=rep(store.p.vec[1],2); py=c(0,1); 
points(px,py,type="l",lwd=1,lty=2); 
elas.sum<-apply(Pelas,c(1,2),sum)*h*h
image(y,y,t(elas.sum),
col=grey((100:500)/500),
xlab="Size year t",ylab="Size year t+1")
contour(y,y,t(elas.sum),method="flattest",add=TRUE,vfont = c("sans serif", "plain"))
#text(locator(1),"c)",cex=1.35,col="white")
#savePlot(filename="e:\\integral\\Onopordum\\FinalFigure3",type="pdf"); 

#plots for B elasticities
win.graph()
par(mfrow=c(2,2),pty="s",bty="l",cex.axis=1.35,cex.lab=1.35)
elas.sum<-apply(BigBelas,4,sum)*h*h	
plot(0:7,elas.sum,type="h",lwd=10,xlab="Age (years)",ylab="Elasticity")
#text(locator(1),"a)",cex=1.35)
elas.sum<-apply(BigBelas,5,sum)*h*h
plot(store.p.vec[1]+cuts.sigma.d,elas.sum,type="h",lwd=2,xlab="Recruit survival intercept",ylab="Elasticity")
#text(locator(1),"b)",cex=1.35)

elas.sum<-apply(BigBelas,c(2,3),sum)*h*h
image(y,y,t(elas.sum),col=grey((100:500)/500),xlab="Size year t",ylab="Size year t+1")
contour(y,y,t(elas.sum),method="flattest",add=TRUE,vfont = c("sans serif", "plain"),nlevels=4)
# text(locator(1),"c)",cex=1.35,col="white")
elas.sum<-apply(BigBelas,1,sum)*h*h
plot(store.p.vec[1]+cuts.sigma.d,elas.sum,type="h",lwd=2,xlab="Adult survival intercept",ylab="Elasticity")
# text(locator(1),"d)",cex=1.35)

#savePlot(filename="e:\\integral\\Onopordum\\FinalFigure4",type="pdf"); 

