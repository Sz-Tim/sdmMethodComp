#============================================================================# 
# R code used to estimate parameters for Onopordum model 
# Note: this file is included to document our fitting procedures. To run 
# the code you would need the data files, which are not provided here.
# Model selection is mostly omitted -- these are just the final fits. 
#
# Beware: this code has proved vulnerable to package updates during our work. 
# Most final calculations for the paper are based on this code using R >= 2.0
# It has been checked as functioning with R 2.1.1 and package versions:
# VR 7.2-20,  glmmML 0.26-3, nlme 3.1-65   
#============================================================================# 
#clear everything, just to be safe 
rm(list=ls(all=TRUE))

require(nlme); require(MASS); require(glmmML); 

#============================================================================# 
#	Fit flowering probability
#============================================================================# 
dataf = data.frame(read.csv("e:\\integral\\onopordum\\big file.fre",header=T))
attach(dataf)

flow<-ifelse(dead==1,NA,0)
flow<-ifelse(is.na(dead),1,flow)
flow<-ifelse(!is.na(lst1),0,flow)

dataf$dead[lst1>0]<-0
dataf$dead[flow==1]<-0
dataf$site= factor(dataf$site)
dataf=data.frame(dataf,flow)
attach(dataf)

fit.flow<-glm(flow~lst+age,family=binomial,data=dataf)
size.fl<-lst[!is.na(lst) & site==1 & flow==1][!is.na(flow)]
size.fl<-size.fl[!is.na(size.fl)]

#============================================================================# 
# Fit seedling size distribution, normal truncated at zero
#============================================================================# 
seedlings<-lst[age==1 & site==1]
seedlings<-seedlings[!is.na(seedlings)]

lik<-function(p){
	lik<-sum(log(dnorm(seedlings,p[1],p[2])/(1-pnorm(0,p[1],p[2]))))
return(-lik)
}
tmp<-optim(c(1,1),lik)
mean.size<-tmp$par[1]
var.size<-tmp$par[2]^2


win.graph(); par(bty="l")
hist(seedlings,col="grey",xlab="Seedling size",main="")

s<-seq(0,7,length=100)
d<-dnorm(s,tmp$par[1],tmp$par[2])/(1-pnorm(0,tmp$par[1],tmp$par[2]))
diff<-s[2]-s[1]
lines(s,d*length(seedlings)/(2*sum(d*diff)))

# overplot normal distribution with same mean and variance 
d<-dnorm(s,mean(seedlings),sd(seedlings))
lines(s,d*length(seedlings)/(2*sum(d*diff)),col="blue")

#============================================================================# 
# Fit growth of surviving plants
#============================================================================# 
# We use GLS to deal with decreasing variance, and test for site effects
# Commented out tests lead to the final model, fitted with REML 
#fit.grow.gls<-gls(lst1~lst+factor(site),na.action=na.omit,weight=varExp(form=~fitted(.)),method="ML")
#fit.grow.gls.0<-gls(lst1~lst+factor(site),na.action=na.omit,method="ML")
#fit.grow.gls.1<-gls(lst1~lst+factor(site),na.action=na.omit,weight=varExp(form=~fitted(.)|site),method="ML")
#anova(fit.grow.gls.0,fit.grow.gls,fit.grow.gls.1)

#test for interaction between site and size
#fit.grow.gls<-gls(lst1~lst+factor(site),na.action=na.omit,weight=varExp(form=~fitted(.)|site),method="ML")
#fit.grow.gls.0<-gls(lst1~lst*factor(site),na.action=na.omit,weight=varExp(form=~fitted(.)|site),method="ML")
#anova(fit.grow.gls.0,fit.grow.gls)

#fit final model by REML
fit.grow.gls<-gls(lst1~lst+factor(site),na.action=na.omit,weight=varExp(form=~fitted(.)|site),data=dataf,verbose=T,
control=glsControl(singular.ok=T))
summary(fit.grow.gls)

var.exp.coef<-as.numeric(fit.grow.gls$modelStruct$varStruct[1])
plot(lst[site==1],lst1[site==1],pch=19)
abline(fit.grow.gls$coef[1],fit.grow.gls$coef[2])

px=seq(from=0,to=9,by=.25); 
mux<-fit.grow.gls$coef[1]+fit.grow.gls$coef[2]*px;
sigmax=(fit.grow.gls$sigma)*exp(var.exp.coef*mux);
points(px,mux+1.96*sigmax,type="l",lty=2); 
points(px,mux-1.96*sigmax,type="l",lty=2); 

mux<-fit.grow.gls$coef[1]+fit.grow.gls$coef[2]*px;
sigmax=(fit.grow.gls$sigma)*exp(var.exp.coef*mux);
points(px,mux+1.96*sigmax,type="l",lty=2); 
points(px,mux-1.96*sigmax,type="l",lty=2); 

#get max and min sizes 
size<-lst[!is.na(lst) & site==1]
minsize<-0.9*min(size); maxsize<-1.1*max(size);
detach(dataf)

#============================================================================# 
# Fit seed production function: negative binomial glm 
#============================================================================# 
data.seeds<-data.frame(read.csv("e:\\integral\\onopordum\\seed production.fre",header=T))
attach(data.seeds)

fit.seeds.nb<-glm.nb(seeds~lst)
summary(fit.seeds.nb)
detach(data.seeds)

#============================================================================# 
# Fit survival: generalized linear mixed model using glmmML
#============================================================================# 
dataf<-data.frame(read.csv("e:\\integral\\onopordum\\deadplno.csv",header=T))
attach(dataf)
surv<-M<0.5
fit.surv<-glmmML(surv~factor(SITE)+LST+AGE,family=binomial,cluster=factor(PLNO))
print.glmmML(fit.surv)
fit.surv.no.var<-glm(surv~factor(SITE)+LST+AGE,family=binomial)

#============================================================================# 
# Store parameter vector
#============================================================================# 
p.vec.names<-rep(NA,15)
p.vec<-rep(0,15);
p.vec[1]<- fit.surv$coef[1]					; p.vec.names[1]<-"Inter survival         ";
p.vec[2]<- fit.surv$coef[3]					; p.vec.names[2]<-"Size slope survival    ";
p.vec[3]<- fit.surv$coef[4]					; p.vec.names[3]<-"Age slope survival     ";
p.vec[4]<- fit.surv$sigma					; p.vec.names[4]<-"sigma survival         ";
p.vec[5]<- fit.flow$coef[1]					; p.vec.names[5]<-"Inter flower           ";
p.vec[6]<- fit.flow$coef[2]					; p.vec.names[6]<-"Size slope flower      ";
p.vec[7]<- fit.flow$coef[3]					; p.vec.names[7]<-"Age slope flower       ";
p.vec[8]<- fit.grow.gls$coef[1]					; p.vec.names[8]<-"ag		          ";
p.vec[9]<- fit.grow.gls$coef[2]					; p.vec.names[9]<-"bg		          ";
p.vec[10]<- fit.grow.gls$sigma^2				; p.vec.names[10]<-"sigma2 growth         ";
p.vec[11]<- fit.seeds.nb$coef[1]  				; p.vec.names[11]<-"intercept seeds       ";
p.vec[12]<- fit.seeds.nb$coef[2]  				; p.vec.names[12]<-"slope seeds           ";
p.vec[13]<- mean.size						; p.vec.names[13]<-"mean kids size        ";
p.vec[14]<- var.size						; p.vec.names[14]<-"sigma2 kids size      ";
p.vec[15]<- var.exp.coef				        ; p.vec.names[15]<-"varExp para           ";

store.p.vec<-p.vec
