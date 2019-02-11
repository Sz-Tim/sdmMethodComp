
sp <- c("barberry", "garlic_mustard")[1]
pkgs <- c("tidyverse", "magrittr", "stringr", "here", "viridis")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
sp_i <- read.csv("data/species_3km.csv") %>% filter(Name==sp)
vs.dir <- paste0("vs/", sp_i$Num, "/")
p <- readRDS(here(vs.dir, "p.rds"))
p$p <- p$germ_x
lam.df <- readRDS(here(vs.dir, "lam_df.rds"))

ipm.diag <- setNames(dir("out/sp1", "Diag_IPM", full.names=T) %>% map(readRDS),
                     str_remove(str_remove(dir("out/sp1", "Diag_IPM"), ".rds"),
                                "Diag_IPM_"))
cad.diag <- setNames(dir("out/sp1", "Diag_CAd", full.names=T) %>% map(readRDS),
                     str_remove(str_remove(dir("out/sp1", "Diag_CAd"), ".rds"),
                                "Diag_CAd_"))
iss.col <- RColorBrewer::brewer.pal(max(length(ipm.diag), length(cad.diag)), "Dark2")
iss.col <- apply(col2rgb(iss.col)/255, 2, function(x) rgb(x[1],x[2],x[3],0.5)) 

z.seq <- seq(p$z.rng[1], p$z.rng[2], length.out=200)
z.mx <- cbind(1, z.seq)
x1 <- seq(min(lam.df$bio10_6), max(lam.df$bio10_6), length.out=200)
x1.mx <- cbind(x1, x1^2)
x2 <- seq(min(lam.df$bio10_prMay), max(lam.df$bio10_prMay), length.out=200)
x2.mx <- cbind(x2, x2^2)



########-------------------------------
## IPM: Size
########
par(mfrow=c(2,3))
plot_sdm_reg(p, ipm.diag, "s_z", 1:length(p$s_z), z.seq, z.mx, iss.col, 
             xlab="Size (t)", ylab="Survival probability", 
             xlim=p$z.rng, ylim=c(0,1))
plot_sdm_reg(p, ipm.diag, "g_z", 1:length(p$g_z), z.seq, z.mx, iss.col, 
             xlab="Size (t)", ylab="Size (t+1)", xlim=p$z.rng, ylim=p$z.rng)
plot_sdm_reg(p, ipm.diag, "fl_z", 1:length(p$fl_z), z.seq, z.mx, iss.col, 
             xlab="Size (t)", ylab="Flowering probability", 
             xlim=p$z.rng, ylim=c(0,1))
plot_sdm_reg(p, ipm.diag, "seed_z", 1:length(p$seed_z), z.seq, z.mx, iss.col, 
             xlab="Size (t)", ylab="Seed production", 
             xlim=p$z.rng, ylim=c(0,6e3))
plot_sdm_reg(p, ipm.diag, "rcr_z", 1:length(p$rcr_z), z.seq, z.mx, iss.col, 
             xlab="Recruit size", ylab="Density", xlim=p$z.rng, ylim=c(0,0.4))
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.col, 1), cex=1.5,
       lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))



########-------------------------------
## IPM: Temperature
########
par(mfrow=c(2,3))
plot_sdm_reg(p, ipm.diag, "s_x", 1:2, x1, x1.mx, iss.col, xlab="Temperature", 
             ylab="Survival probability", xlim=range(x1), ylim=c(0,1))
plot_sdm_reg(p, ipm.diag, "g_x", 1:2, x1, x1.mx, iss.col, xlab="Temperature", 
             ylab="Size (t+1)", xlim=range(x1), ylim=c(-8,1))
plot_sdm_reg(p, ipm.diag, "fl_x", 1:2, x1, x1.mx, iss.col, xlab="Temperature", 
             ylab="Flowering probability", xlim=range(x1), ylim=c(0,1))
plot_sdm_reg(p, ipm.diag, "seed_x", 1:2, x1, x1.mx, iss.col, xlab="Temperature", 
             ylab="Seed production", xlim=range(x1), ylim=c(-2,2))
plot_sdm_reg(p, ipm.diag, "germ_x", 1:3, x1, cbind(1, x1.mx), iss.col, 
             xlab="Temperature", ylab="Germination probability", 
             xlim=range(x1), ylim=c(0,1))
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.col, 1), cex=1.5,
       lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))



########-------------------------------
## IPM: Precipitation
########
par(mfrow=c(2,3))
plot_sdm_reg(p, ipm.diag, "s_x", 3:4, x2, x2.mx, iss.col, xlab="Precipitation", 
             ylab="Survival probability", xlim=range(x2), ylim=c(0,1))
plot_sdm_reg(p, ipm.diag, "g_x", 3:4, x2, x2.mx, iss.col, xlab="Precipitation", 
             ylab="Size (t+1)", xlim=range(x2), ylim=c(-8,1))
plot_sdm_reg(p, ipm.diag, "fl_x", 3:4, x2, x2.mx, iss.col, xlab="Precipitation", 
             ylab="Flowering probability", xlim=range(x2), ylim=c(0,1))
plot_sdm_reg(p, ipm.diag, "seed_x", 3:4, x2, x2.mx, iss.col, xlab="Precipitation", 
             ylab="Seed production", xlim=range(x2), ylim=c(-2,2))
plot_sdm_reg(p, ipm.diag, "germ_x", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
             xlab="Precipitation", ylab="Germination probability", 
             xlim=range(x2), ylim=c(0,1))
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.col, 1), cex=1.5,
       lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))



########-------------------------------
## CAd: Temperature
########
par(mfrow=c(2,3))
plot_sdm_reg(p, cad.diag, "s.M", 1:3, x1, cbind(1, x1.mx), iss.col, 
             xlab="Temperature", ylab="Juvenile survival", 
             xlim=range(x1), ylim=c(0,1))
# plot_sdm_reg(p, cad.diag, "s.N", 1:3, x1, cbind(1, x1.mx), iss.col,
#              xlab="Temperature", ylab="Adult survival",
#              xlim=range(x1), ylim=c(0,1))
plot_sdm_reg(p, cad.diag, "K", 1:3, x1, cbind(1, x1.mx), iss.col, 
             xlab="Temperature", ylab="Carrying capacity", 
             xlim=range(x1), ylim=c(0,3e3))
plot_sdm_reg(p, cad.diag, "p.f", 1:3, x1, cbind(1, x1.mx), iss.col, 
             xlab="Temperature", ylab="Flowering probability", 
             xlim=range(x1), ylim=c(0,1))
plot_sdm_reg(p, cad.diag, "mu", 1:3, x1, cbind(1, x1.mx), iss.col, 
             xlab="Temperature", ylab="Seed production", 
             xlim=range(x1), ylim=c(0,6e3))
plot_sdm_reg(p, cad.diag, "p", 1:3, x1, cbind(1, x1.mx), iss.col, 
             xlab="Temperature", ylab="Establishment probability", 
             xlim=range(x1), ylim=c(0,1))
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=rep(1, length(iss.col)), col=iss.col, cex=1.2,
       lwd=rep(3, length(iss.col)), legend=names(cad.diag))



########-------------------------------
## CAd: Precipitation
########
par(mfrow=c(2,3))
plot_sdm_reg(p, cad.diag, "s.M", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
             xlab="Precipitation", ylab="Juvenile survival", 
             xlim=range(x2), ylim=c(0,1))
# plot_sdm_reg(p, cad.diag, "s.N", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
#              xlab="Precipitation", ylab="Adult survival", 
#              xlim=range(x2), ylim=c(0,1))
plot_sdm_reg(p, cad.diag, "K", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
             xlab="Precipitation", ylab="Carrying capacity", 
             xlim=range(x2), ylim=c(0,3e3))
plot_sdm_reg(p, cad.diag, "p.f", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
             xlab="Precipitation", ylab="Flowering probability", 
             xlim=range(x2), ylim=c(0,1))
plot_sdm_reg(p, cad.diag, "mu", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
             xlab="Precipitation", ylab="Seed production", 
             xlim=range(x2), ylim=c(0,6e3))
plot_sdm_reg(p, cad.diag, "p", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
             xlab="Precipitation", ylab="Establishment probability", 
             xlim=range(x2), ylim=c(0,1))
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=rep(1, length(iss.col)), col=iss.col, cex=1.2,
       lwd=rep(3, length(iss.col)), legend=names(cad.diag))







