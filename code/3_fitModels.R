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
sampling.issue <- c("none", "error", "geog", "bias")[1]
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
n_sim <- 15  # number of simulations per sample (mechanistic only)
v <- c("(Intercept)"=0, "size"=0, "size2"=0, "size3"=0, 
       "temp"=0, "temp2"=0, "prec"=0, "prec2"=0, 
       "pOpn"=0, "pOth"=0, "pDec"=0, "pEvg"=0, "pMxd"=0)
m.full <- paste(names(v)[-1], collapse=" + ")
n_z <- sum(grepl("size", names(v))) + 1  # include intercept
n_x <- length(v) - n_z
n_z <- rep(list(n_z), 4)
names(n_z) <- c("s", "g", "fl", "seed")
n_x <- rep(list(n_x), 4)
names(n_x) <- c("s", "g", "fl", "seed")
X <- list(s=as.matrix(env.in[,1:n_x$s]), g=as.matrix(env.in[,1:n_x$g]),
          fl=as.matrix(env.in[,1:n_x$fl]), seed=as.matrix(env.in[,1:n_x$seed]))


########
## fit models
########
# MaxEnt
# CA
# IPM
S.f <- U.f <- vector("list", length(O_IPM))
for(i in 1:length(O_IPM)) {
  # separate data from MuMIn::dredge()
  O_IPM.i.s <- filter(O_IPM[[i]], !is.na(surv))
  O_IPM.i.g <- filter(O_IPM[[i]], !is.na(sizeNext) & !is.na(size))
  O_IPM.i.fl <- filter(O_IPM[[i]], !is.na(fl))
  O_IPM.i.seed <- filter(O_IPM[[i]], !is.na(seed))
  
  # global models
  options(na.action="na.fail")
  s.m <- glm(as.formula(paste("surv ~", m.full, collapse="")), 
             data=O_IPM.i.s, family="binomial")
  g.m <- lm(as.formula(paste("sizeNext ~", m.full, collapse="")), 
            data=O_IPM.i.g)
  fl.m <- glm(as.formula(paste("fl ~", m.full, collapse="")), 
              data=O_IPM.i.fl, family="binomial")
  seed.m <- glm(as.formula(paste("seed ~", m.full, collapse="")), 
                data=O_IPM.i.seed, family="poisson")
  
  # optimal models
  s.opt <- get.models(dredge(s.m), subset=1)[[1]]
  g.opt <- get.models(dredge(g.m), subset=1)[[1]]
  fl.opt <- get.models(dredge(fl.m), subset=1)[[1]]
  seed.opt <- get.models(dredge(seed.m), subset=1)[[1]]
  
  # store coefficients
  vars.opt <- list(s=coef(s.opt),
                   g=coef(g.opt),
                   fl=coef(fl.opt),
                   seed=coef(seed.opt))
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
  p.ipm$g_sig <- summary(g.opt)$sigma
  p.ipm$fl_z <- vars.ls$fl[1:n_z$fl]
  p.ipm$fl_x <- vars.ls$fl[(n_z$fl+1):length(v)]
  p.ipm$seed_z <- vars.ls$seed[1:n_z$seed]
  p.ipm$seed_x <- vars.ls$seed[(n_z$seed+1):length(v)]
  p.ipm$rcr_z <- filter(O_IPM[[i]], is.na(size)) %>% 
    summarise(mn=mean(sizeNext), sd=sd(sizeNext)) %>% unlist
  
  U.f[[i]] <- fill_IPM_matrices(n.cell, buffer=0.75, discrete=1, p.ipm, 
                                n_z, n_x, X, sdd.pr, env.in$id, verbose=T)
  
  #---- Realization: use estimated slopes to generate simulated data
  S.i <- vector("list", n_sim)
  for(s in 1:n_sim) {
    S.i[[s]] <- simulate_data(n.cell, U$lo, U$hi, p.ipm, X, n_z, 
                              sdd.pr, U$sdd.j, verbose=T)
    cat("||-- Finished simulation", s, "\n||--\n")
  }
  S.f[[i]] <- summarize_IPM_simulations(S.i, p.ipm$tmax)
}







