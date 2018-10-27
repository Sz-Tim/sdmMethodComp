### R code for Figure 5 of Ellner and Rees, Am Nat (2005 or thereabouts). 
### The code can be copied and pasted directly into R, or sourced.

### This script uses functions for the Onopordum model defined in RunOnopordumModel.R, 
### so you need to source that file first or this one won't run.  

### The calculations are just nested loops in which the model is
### iterated to convergence for different values of the parameters,  
### to find the stable growth rate or the equilibrium population size.  

### So that it runs quickly this code uses very few mesh points for size and quality
### For accurate results we recommend using 
	# n.big.matrix <- 50; n.sigma <- 40; 
### You may need add --max-mem-size=2G to the end of the target line of your R shortcut 
### (in MS Windows) to use this number of mesh points.

#============================================================================================#
# Parameter values  
#============================================================================================#
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
store.p.vec=p.vec; base.p.vec=p.vec;  
tol=1e-6; 

############### Set up matrices for iteration 
#P=survival-growth, B=fecundity

P<-array(NA,dim=c(n.big.matrix,n.big.matrix,n.age,n.sigma))
B<-array(NA,dim=c(n.big.matrix,n.big.matrix,n.age,n.sigma))

L<-minsize; U<-maxsize;
n<-n.big.matrix
	
# boundary points b and mesh points y
b<-L+c(0:n)*(U-L)/n;
y<-0.5*(b[1:n]+b[2:(n+1)]);

h<-y[2]-y[1]

prob.offspring<-dnorm(y,p.vec[13],sqrt(p.vec[14]))
prob.offspring<-h*prob.offspring/sum(prob.offspring*h)

#points for numerical integration of survival intercepts
cuts.sigma.d<-seq(-4*p.vec[4],4*p.vec[4],length=n.sigma);
prob.sigma.d<-dnorm(cuts.sigma.d)/sum(dnorm(cuts.sigma.d));

# Create P and B matrices 
for (inter in 1:n.sigma){
 for (age in 1:n.age){
	p.vec[1]<-store.p.vec[1]+cuts.sigma.d[inter]
	P[,,age,inter]<-(U-L)*t(outer(y,y,pxy,age=age,params=p.vec))/n
	B[,,age,inter]<-(U-L)*t(outer(y,y,fxy,age=age,params=p.vec))/n
}
}

# Create B.short 
B.short<-apply(B,c(2,3,4),sum)
dist.sigma.d.by.offspring.size<-outer(prob.offspring,prob.sigma.d)	

############# vary sigma_s and calculate lambda and equilibrium pop size. 
n.sigmas=25
sigma.s<-seq(0.001,2,length=n.sigmas)
lams<-rep(NA,n.sigmas)
equ.pop.size<-rep(NA,n.sigmas)

for(i in 1:n.sigmas){
p.vec[4]<-sigma.s[i];

#points for numerical integration of survival intercepts
cuts.sigma.d<-seq(-4*p.vec[4],4*p.vec[4],length=n.sigma);
prob.sigma.d<-dnorm(cuts.sigma.d)/sum(dnorm(cuts.sigma.d));

# Create P and B matrices 
for (inter in 1:n.sigma){
 for (age in 1:n.age){
	p.vec[1]<-store.p.vec[1]+cuts.sigma.d[inter]
	P[,,age,inter]<-(U-L)*t(outer(y,y,pxy,age=age,params=p.vec))/n
	B[,,age,inter]<-(U-L)*t(outer(y,y,fxy,age=age,params=p.vec))/n
}
}

# Create B.short and BigB 
B.short<-apply(B,c(2,3,4),sum)
dist.sigma.d.by.offspring.size<-outer(prob.offspring,prob.sigma.d)	

qmax=1; lam=1; iter=1;
Nt<-array(1,dim=c(n.big.matrix,n.age,n.sigma));
Nt=Nt/sum(Nt); 

while (qmax>tol){
	Nt1=iteration(Nt,use.B.short=T)
	qmax=sum(abs(Nt1-lam*Nt)); 
	lam=sum(Nt1); 
	#cat(iter,sum(Nt),"  ",sum(Nt1),"  ",lam, qmax,"\n")
	iter=iter+1; 
	Nt = Nt1/lam;
}
lams[i]<-lam

if(lams[i]>1){
	Nt<-array(1,dim=c(n.big.matrix,n.age,n.sigma));
	Nt=Nt/sum(Nt); 
	iter=1;
	while (iter<100){
		Nt1=iteration(Nt,use.B.short=T,DI=F)  
		#cat(iter,sum(Nt),"  ",sum(Nt1),"  ", qmax,"\n")
		iter=iter+1; 
		Nt = Nt1;
	}
	equ.pop.size[i]<-sum(Nt)

} else equ.pop.size[i]<-0
cat(i,"  ",sigma.s[i],"  ",lam,"  ",equ.pop.size[i],"\n")
} 

###### Vary the mean offspring quality, for 3 possible values of quality variance 
sigma.inters=base.p.vec[4]*c(1,0.25,4);  
mu.inters=seq(from=base.p.vec[1]-2*base.p.vec[4],to=base.p.vec[1]+2*base.p.vec[4],length=25); 

lamvals=matrix(0,25,3) 
ntvals=lamvals; 

for(i in 1:25) {
for(j in 1:3) {
store.p.vec=base.p.vec;
store.p.vec[1]=mu.inters[i];
store.p.vec[4]=sigma.inters[j];  
p.vec=store.p.vec; 
	
prob.offspring<-dnorm(y,p.vec[13],sqrt(p.vec[14]))
prob.offspring<-h*prob.offspring/sum(prob.offspring*h)

#points for numerical integration of survival intercepts
cuts.sigma.d<-seq(-4*p.vec[4],4*p.vec[4],length=n.sigma);
prob.sigma.d<-dnorm(cuts.sigma.d)/sum(dnorm(cuts.sigma.d));

# Create P and B matrices 
for (inter in 1:n.sigma){
 for (age in 1:n.age){
	p.vec[1]<-store.p.vec[1]+cuts.sigma.d[inter]
	P[,,age,inter]<-(U-L)*t(outer(y,y,pxy,age=age,params=p.vec))/n
	B[,,age,inter]<-(U-L)*t(outer(y,y,fxy,age=age,params=p.vec))/n
}
}

# Create B.short 
B.short<-apply(B,c(2,3,4),sum)
dist.sigma.d.by.offspring.size<-outer(prob.offspring,prob.sigma.d)	

################ iterate model using B.short to get lam and stable.dist 
Nt<-array(1,dim=c(n.big.matrix,n.age,n.sigma));
Nt=Nt/sum(Nt); 

qmax=1; lam=1; iter=1;

while (qmax>tol){
	Nt1=iteration(Nt,use.B.short=T)
	qmax=sum(abs(Nt1-lam*Nt)); 
	lam=sum(Nt1); 
	iter=iter+1; 
	Nt = Nt1/lam;
}
lamvals[i,j]=lam; 

if(lamvals[i,j]>1) {

	Nt<-array(1,dim=c(n.big.matrix,n.age,n.sigma));
	Nt=Nt/sum(Nt); 
	iter=1; while (iter<100){
		Nt1=iteration(Nt,use.B.short=T,DI=F)  
		iter=iter+1; 
		Nt = Nt1;
	}
	ntvals[i,j]<-sum(Nt)

} else ntvals[i,j]=0

cat(i,j, lamvals[i,j], ntvals[i,j],"\n")

}}

win.graph()
par(mfrow=c(2,2),bty="l",pty="s", mar=c(3,6,3,3),
cex.axis=1.5,cex.lab=1.5)
plot(sigma.s,lams,type="l",xlim=c(0,2),ylim=c(0.98,1.04),lwd=2,xlab="Std Dev Offspring Quality",ylab=expression("Finite rate of increase, "*lambda))
abline(v=store.p.vec[4],lty=2)
points(0,1.035838,pch=19,cex=1.75)
text(0.25,1.04,"a)",cex=1.35) 
ff=approxfun(sigma.s,lams); 
xx=base.p.vec[4]; yy=ff(xx);  
points(xx,yy,pch=2,cex=2,lwd=2); 

plot(sigma.s,equ.pop.size/40,xlim=c(0,2),ylim=c(0,5), type="l",lwd=2,xlab="Std Dev Offspring Quality",
ylab=expression("Equilibrium density, plants  " *m^-2))
abline(v=store.p.vec[4],lty=2)
points(0,189.7001/40,pch=19,cex=1.75)
text(0.25,5,"b)",cex=1.35) 
ff=approxfun(sigma.s,equ.pop.size/40); 
xx=base.p.vec[4]; yy=ff(xx);  
points(xx,yy,pch=2,cex=2,lwd=2); 

plot(mu.inters,lamvals[,1],type="l",col="black",xlab="Mean Offspring Quality",lwd=2,
ylab=expression("Finite rate of increase, "*lambda))
points(mu.inters[13],lamvals[13,1],pch=2,cex=1.7,lwd=2); 
text(-2.5,1.2,"c)",cex=1.35) 

plot(mu.inters,ntvals[,1]/40,type="l",col="black",xlab="Mean Offspring Quality",lwd=2,
ylab=expression("Equilibrium density, plants  " *m^-2)) 
points(mu.inters[13],ntvals[13,1]/40,pch=2,cex=1.7,lwd=2); 
text(-2.5,7,"d)",cex=1.35) 

