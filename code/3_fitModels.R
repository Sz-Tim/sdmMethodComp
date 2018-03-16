# 3: Fit models from samples
# Comparison of SDM approaches
# Tim Szewczyk

# This script uses the samples generated in `2_sampleDistribution.R` to fit a
# correlative MaxEnt SDM, a mechanistic CA-style SDM, and a mechanistic 
# IPM-style SDM. Model selection is performed for each model

########
## Setup
########
# file specifications
sp <- "sp1"
sampling.issue <- c("none", "error", "geog", "bias")[4]
modeling.issue <- c("none")[1]

# load workspace
pkgs <- c("gbPopMod", "tidyverse", "magrittr", "MuMIn")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
source("code/fn_IPM.R"); source("code/fn_aux.R"); source("code/fn_sim.R")
p <- readRDS(paste0("out/", sp, "_p.rds"))
S <- readRDS(paste0("out/", sp, "_S.rds"))
U <- readRDS(paste0("out/", sp, "_U.rds"))
sdd.pr <- readRDS(paste0("out/", sp, "_sdd.rds"))
lam.df <- readRDS(paste0("out/", sp, "_lam_df.rds"))
env.in <- readRDS(paste0("out/", sp, "_env_in.rds"))
n.cell <- nrow(env.in)
O_Mx <- readRDS(paste0("out/", sp, "_O_Mx_", sampling.issue, ".rds"))
O_CA <- readRDS(paste0("out/", sp, "_O_CA_", sampling.issue, ".rds"))
O_IPM <- readRDS(paste0("out/", sp, "_O_IPM_", sampling.issue, ".rds"))


########
## Set model details
########
n_sim <- 10  # number of simulations per sample (mechanistic only)
v <- c("(Intercept)"=0, "size"=0, "size2"=0,# "size3"=0, 
       "temp"=0, "temp2"=0, "prec"=0, "prec2"=0, 
       "pOpn"=0, "pOth"=0, "pDec"=0)#, "pEvg"=0, "pMxd"=0)
m.full <- paste(names(v)[-1], collapse=" + ")
n_z <- rep(list(sum(grepl("size", names(v))) + 1), 4)# include intercept
n_x <- rep(list(length(v) - n_z[[1]]), 4)
names(n_z) <- names(n_x) <- c("s", "g", "fl", "seed")
X <- map(n_x, ~as.matrix(env.in[,1:.]))

########
## fit models
########
##--- MaxEnt
##--- CA
##--- IPM
S.f <- U.f <- vector("list", length(O_IPM))
for(i in 1:length(O_IPM)) {
  # separate data for MuMIn::dredge()
  O_IPM.i.s <- filter(O_IPM[[i]], !is.na(surv))
  O_IPM.i.g <- filter(O_IPM[[i]], !is.na(sizeNext) & !is.na(size))
  O_IPM.i.fl <- filter(O_IPM[[i]], !is.na(fl))
  O_IPM.i.seed <- filter(O_IPM[[i]], !is.na(seed))
  
  # global models
  options(na.action="na.fail")
  global.m <- list(s=glm(as.formula(paste("surv ~", m.full, collapse="")), 
                         data=O_IPM.i.s, family="binomial"),
                   g=lm(as.formula(paste("sizeNext ~", m.full, collapse="")), 
                        data=O_IPM.i.g),
                   fl=glm(as.formula(paste("fl ~", m.full, collapse="")), 
                          data=O_IPM.i.fl, family="binomial"),
                   seed=glm(as.formula(paste("seed ~", m.full, collapse="")), 
                            data=O_IPM.i.seed, family="poisson"))
  
  # optimal models
  opt.m <- map(global.m, ~get.models(dredge(.), subset=1)[[1]])
  
  # store coefficients
  vars.opt <- map(opt.m, coef)
  vars.ls <- list(s=v, g=v, fl=v, seed=v)
  vars.ls$s[names(vars.opt$s)] <- vars.opt$s
  vars.ls$g[names(vars.opt$g)] <- vars.opt$g
  vars.ls$fl[names(vars.opt$fl)] <- vars.opt$fl
  vars.ls$seed[names(vars.opt$seed)] <- vars.opt$seed
  
  # update parameters
  p.ipm <- p
  p.ipm$s_z <- vars.ls$s[1:n_z$s]
  p.ipm$s_x <- vars.ls$s[(n_z$s+1):length(v)]
  p.ipm$g_z <- vars.ls$g[1:n_z$g]
  p.ipm$g_x <- vars.ls$g[(n_z$g+1):length(v)]
  p.ipm$g_sig <- summary(opt.m$g)$sigma
  p.ipm$fl_z <- vars.ls$fl[1:n_z$fl]
  p.ipm$fl_x <- vars.ls$fl[(n_z$fl+1):length(v)]
  p.ipm$seed_z <- vars.ls$seed[1:n_z$seed]
  p.ipm$seed_x <- vars.ls$seed[(n_z$seed+1):length(v)]
  p.ipm$rcr_z <- filter(O_IPM[[i]], is.na(size)) %>% 
    summarise(mn=mean(sizeNext), sd=sd(sizeNext)) %>% unlist
  
  # use estimated slopes to fill IPM matrix
  U.f[[i]] <- fill_IPM_matrices(n.cell, buffer=0.75, discrete=1, p.ipm, 
                                n_z, n_x, X, sdd.pr, env.in$id)
  
  # use estimated slopes to generate simulated data
  Si <- vector("list", n_sim)
  cat("||-- Starting simulations\n")
  for(s in 1:n_sim) {
    Si[[s]] <- simulate_data(n.cell, U$lo, U$hi, p.ipm, X, n_z, sdd.pr, U$sdd.j)
    cat("||-- Finished simulation", s, "\n")
  }
  S.f[[i]] <- summarize_IPM_simulations(Si, p.ipm$tmax)
  rm(Si)
  cat("\n  Finished dataset", i, "\n")
}
out <- summarize_IPM_samples(U.f, S.f)

P_IPM <- lam.df %>% select("x", "y", "x_y", "lat", "lon", "id", "id.inbd") %>% 
  mutate(lambda.f=apply(out$Uf$IPM.mn, 3, function(x) Re(eigen(x)$values[1])),
         lam.S.f=rowMeans(out$Sf$N_sim.mn[,(-3:0)+p.ipm$tmax]/
                            (out$Sf$N_sim.mn[,(-4:-1)+p.ipm$tmax]+0.0001)),
         nSeed.f=out$Sf$nSd.mn[,p.ipm$tmax], 
         D.f=out$Sf$D.mn[,p.ipm$tmax],
         B0.f=out$Sf$B.mn[,1], 
         Btmax.f=out$Sf$B.mn[,p.ipm$tmax+1],
         N.S.f=out$Sf$N_tot.mn, 
         Surv.S.f=out$Sf$N_surv.mn, 
         Rcr.S.f=out$Sf$N_rcr.mn,
         nSdStay.f=nSeed.f*(1-p.ipm$p_emig), 
         nSdLeave.f=nSeed.f*p.ipm$p_emig,
         N.U.f=apply(out$Uf$Nt.mn[,,p.ipm$tmax],2,sum), 
         lam.U.f=out$Uf$lam.mn[,p.ipm$tmax-1])


saveRDS(P_IPM, paste0("out/", sp, "_P_IPM_", sampling.issue, "_",
                      modeling.issue, ".rds"))





