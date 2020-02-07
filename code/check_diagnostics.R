
sp <- c("barberry", "garlic_mustard")[2]
spType <- c("shrub", "biennial")[(sp=="garlic_mustard")+1]
write.plots <- TRUE
by.issue <- FALSE
pkgs <- c("tidyverse", "magrittr", "stringr", "here", "viridis", "gridExtra")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
sp_i <- read.csv("data/species_3km.csv") %>% filter(Name==sp)
vs.dir <- paste0("vs/", sp_i$Num, "/")
out.dir <- paste0("out/", sp_i$Num, "/")
p <- readRDS(here(vs.dir, "p.rds"))
p$g.D <- p$germ_x
lam.df <- readRDS(here(vs.dir, "lam_df.rds"))
issues <- c("none", "noise", "nonEq", "sampBias", "wrongCov")

ipm.diag <- setNames(dir(out.dir, "Diag_IPM", full.names=T) %>% map(readRDS),
                     str_remove(str_remove(dir(paste0("out/", sp_i$Num), 
                                               "Diag_IPM"), ".rds"), "Diag_IPM_"))
ipm.diag <- ipm.diag[issues]
cad.diag <- setNames(dir(out.dir, "Diag_CAd", full.names=T) %>% map(readRDS),
                     str_remove(str_remove(dir(paste0("out/", sp_i$Num), 
                                               "Diag_CAd"), ".rds"), "Diag_CAd_"))
cad.diag <- cad.diag[issues]
mxe.auc <- read_csv(paste0(out.dir, "out_AUC.csv"))
mxe.roc <- read_csv(paste0(out.dir, "out_ROC.csv"))


iss.COL <- RColorBrewer::brewer.pal(max(length(ipm.diag), 
                                             length(cad.diag)), "Dark2")
iss.col <- apply(col2rgb(iss.COL)/255, 2, function(x) rgb(x[1],x[2],x[3],0.1)) 

z.seq <- seq(p$z.rng[1], p$z.rng[2], length.out=200)
z.mx <- cbind(1, z.seq)
x1 <- seq(min(lam.df$bio10_6), max(lam.df$bio10_6), length.out=200)
x1.mx <- cbind(x1, x1^2)
x2 <- seq(min(lam.df$bio10_prMay), max(lam.df$bio10_prMay), length.out=200)
x2.mx <- cbind(x2, x2^2)



########-------------------------------
## IPM: Slope comparisons
########

param.ls <- list(c("s_z", "s_x"), c("g_z", "g_x"), 
                 c("fl_z", "fl_x"), c("seed_z", "seed_x"), 
                 c("germ_x"))
ipm.slopes.ls <- ipm.true.ls <- gg.p <- vector("list", length(param.ls))
reg_vars <- c("(Intercept)"="(Intercept)", "size"="size",
              "minTemp"="bio10_6", "minTemp^2"="bio10_6_sq", 
              "precMay"="bio10_prMay", "precMay^2"="bio10_prMay_sq")
paramTypes <- c("(a) survival", "(b) growth", "(c) flowering", 
                "(d) seeds", "(e) germination")
for(i in seq_along(param.ls)) {
  params <- param.ls[[i]]
  params <- setNames(params, str_split_fixed(params, "_", 2)[,2])
  paramType <- paramTypes[i]
  ipm.slopes <- map_dfr(params, ~plot_sdm_slopes(p, ipm.diag, .)$slope.df) %>%
    mutate(variable=factor(variable, levels=reg_vars, labels=names(reg_vars)),
           param=paramType,
           scenario=factor(scenario, 
                           levels=c("none", "noise", "sampBias", "nonEq"),
                           labels=c("Ideal", "Measurement error",
                                    "Sampling bias", "Non-equilibrium")))
  ipm.true <- map_dfr(params, ~plot_sdm_slopes(p, ipm.diag, .)$true.df) %>%
    mutate(variable=factor(variable, levels=reg_vars, labels=names(reg_vars)),
           param=paramType)
  ipm.slopes.ls[[i]] <- ipm.slopes
  ipm.true.ls[[i]] <- ipm.true 
  gg.p[[i]] <- ggplot() + 
    scale_colour_manual(values=iss.COL) + 
    scale_fill_manual(values=iss.col) + 
    geom_density(data=ipm.slopes, aes(x=value, colour=scenario, fill=scenario)) +
    geom_vline(data=ipm.true, aes(xintercept=value), linetype=2) +
    facet_wrap(~variable, scales="free", drop=F, nrow=1, labeller=label_parsed) +
    theme(panel.grid=element_blank(),
          axis.text=element_text(size=5),
          axis.title=element_text(size=8)) + 
    ggtitle(paste0(paramType, ": ", spType))
}
gg.grid <- grid.arrange(grobs=gg.p, ncol=1)
ggsave(paste0("figs/diag/slope_dens_", spType, ".pdf"), gg.grid, width=9, height=10)






########-------------------------------
## IPM: Size
########
if(by.issue) {
  for(i in 1:length(issues)) {
    f.name <- paste0("figs/diag/byIssue/", sp_i$Num, "_IPM_Size_", issues[i], ".jpg")
    if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
    par(mfrow=c(2,3))
    plot_sdm_reg(p, ipm.diag[i], "s_z", 1:length(p$s_z), z.seq, z.mx, iss.col[i], 
                 xlab="Size (t)", ylab="Survival probability", 
                 xlim=p$z.rng, ylim=c(0,1))
    plot_sdm_reg(p, ipm.diag[i], "g_z", 1:length(p$g_z), z.seq, z.mx, iss.col[i], 
                 xlab="Size (t)", ylab="Size (t+1)", xlim=p$z.rng, ylim=p$z.rng)
    plot_sdm_reg(p, ipm.diag[i], "fl_z", 1:length(p$fl_z), z.seq, z.mx, iss.col[i], 
                 xlab="Size (t)", ylab="Flowering probability", 
                 xlim=p$z.rng, ylim=c(0,1))
    plot_sdm_reg(p, ipm.diag[i], "seed_z", 1:length(p$seed_z), z.seq, z.mx, iss.col[i], 
                 xlab="Size (t)", ylab="Seed production", 
                 xlim=p$z.rng, ylim=c(0,1e3))
    plot_sdm_reg(p, ipm.diag[i], "rcr_z", 1:length(p$rcr_z), z.seq, z.mx, iss.col[i], 
                 xlab="Recruit size", ylab="Density", xlim=p$z.rng, ylim=c(0,0.5))
    plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
    legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.COL, 1), cex=1.5,
           lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))
    if(write.plots) dev.off()
  }
} else {
  f.name <- paste0("figs/diag/", sp_i$Num, "_IPM_Size.jpg")
  if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
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
               xlab="Recruit size", ylab="Density", xlim=p$z.rng, ylim=c(0,0.5))
  plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
  legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.COL, 1), cex=1.5,
         lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))
  if(write.plots) dev.off()
}



########-------------------------------
## IPM: Temperature
########
if(by.issue) {
  for(i in 1:length(issues)) {
    f.name <- paste0("figs/diag/byIssue/", sp_i$Num, "_IPM_Temp_", issues[i], ".jpg")
    if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
    par(mfrow=c(2,3))
    plot_sdm_reg(p, ipm.diag[i], "s_x", 1:2, x1, x1.mx, iss.col[i], xlab="Min Temp: bio6", 
                 ylab="Survival probability", xlim=range(x1), ylim=c(0,1))
    plot_sdm_reg(p, ipm.diag[i], "g_x", 1:2, x1, x1.mx, iss.col[i], xlab="Min Temp: bio6", 
                 ylab="Size (t+1)", xlim=range(x1), ylim=c(-8,4))
    plot_sdm_reg(p, ipm.diag[i], "fl_x", 1:2, x1, x1.mx, iss.col[i], xlab="Min Temp: bio6", 
                 ylab="Flowering probability", xlim=range(x1), ylim=c(0,1))
    plot_sdm_reg(p, ipm.diag[i], "seed_x", 1:2, x1, x1.mx, iss.col[i], xlab="Min Temp: bio6", 
                 ylab="Seed production", xlim=range(x1), ylim=c(-2,2))
    plot_sdm_reg(p, ipm.diag[i], "germ_x", 1:3, x1, cbind(1, x1.mx), iss.col[i], 
                 xlab="Min Temp: bio6", ylab="Germination probability", 
                 xlim=range(x1), ylim=c(0,1))
    plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
    legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.COL, 1), cex=1.5,
           lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))
    if(write.plots) dev.off()
  }
} else {
  f.name <- paste0("figs/diag/", sp_i$Num, "_IPM_Temp.jpg")
  if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
  par(mfrow=c(2,3))
  plot_sdm_reg(p, ipm.diag, "s_x", 1:2, x1, x1.mx, iss.col, xlab="Min Temp: bio6", 
               ylab="Survival probability", xlim=range(x1), ylim=c(0,1))
  plot_sdm_reg(p, ipm.diag, "g_x", 1:2, x1, x1.mx, iss.col, xlab="Min Temp: bio6", 
               ylab="Size (t+1)", xlim=range(x1), ylim=c(-8,4))
  plot_sdm_reg(p, ipm.diag, "fl_x", 1:2, x1, x1.mx, iss.col, xlab="Min Temp: bio6", 
               ylab="Flowering probability", xlim=range(x1), ylim=c(0,1))
  plot_sdm_reg(p, ipm.diag, "seed_x", 1:2, x1, x1.mx, iss.col, xlab="Min Temp: bio6", 
               ylab="Seed production", xlim=range(x1), ylim=c(-2,2))
  plot_sdm_reg(p, ipm.diag, "germ_x", 1:3, x1, cbind(1, x1.mx), iss.col, 
               xlab="Min Temp: bio6", ylab="Germination probability", 
               xlim=range(x1), ylim=c(0,1))
  plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
  legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.COL, 1), cex=1.5,
         lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))
  if(write.plots) dev.off()
}



########-------------------------------
## IPM: Precipitation
########
if(by.issue) {
  for(i in 1:length(issues)) {
    f.name <- paste0("figs/diag/byIssue/", sp_i$Num, "_IPM_Precip_", issues[i], ".jpg")
    if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
    par(mfrow=c(2,3))
    plot_sdm_reg(p, ipm.diag[i], "s_x", 3:4, x2, x2.mx, iss.col[i], xlab="May Precip", 
                 ylab="Survival probability", xlim=range(x2), ylim=c(0,1))
    plot_sdm_reg(p, ipm.diag[i], "g_x", 3:4, x2, x2.mx, iss.col[i], xlab="May Precip", 
                 ylab="Size (t+1)", xlim=range(x2), ylim=c(-8,2))
    plot_sdm_reg(p, ipm.diag[i], "fl_x", 3:4, x2, x2.mx, iss.col[i], xlab="May Precip", 
                 ylab="Flowering probability", xlim=range(x2), ylim=c(0,1))
    plot_sdm_reg(p, ipm.diag[i], "seed_x", 3:4, x2, x2.mx, iss.col[i], xlab="May Precip", 
                 ylab="Seed production", xlim=range(x2), ylim=c(-2,2))
    plot_sdm_reg(p, ipm.diag[i], "germ_x", c(1,4,5), x2, cbind(1, x2.mx), iss.col[i], 
                 xlab="May Precip", ylab="Germination probability", 
                 xlim=range(x2), ylim=c(0,1))
    plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
    legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.COL, 1), cex=1.5,
           lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))
    if(write.plots) dev.off()
  }
} else {
  f.name <- paste0("figs/diag/", sp_i$Num, "_IPM_Precip.jpg")
  if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
  par(mfrow=c(2,3))
  plot_sdm_reg(p, ipm.diag, "s_x", 3:4, x2, x2.mx, iss.col, xlab="May Precip", 
               ylab="Survival probability", xlim=range(x2), ylim=c(0,1))
  plot_sdm_reg(p, ipm.diag, "g_x", 3:4, x2, x2.mx, iss.col, xlab="May Precip", 
               ylab="Size (t+1)", xlim=range(x2), ylim=c(-8,2))
  plot_sdm_reg(p, ipm.diag, "fl_x", 3:4, x2, x2.mx, iss.col, xlab="May Precip", 
               ylab="Flowering probability", xlim=range(x2), ylim=c(0,1))
  plot_sdm_reg(p, ipm.diag, "seed_x", 3:4, x2, x2.mx, iss.col, xlab="May Precip", 
               ylab="Seed production", xlim=range(x2), ylim=c(-2,2))
  plot_sdm_reg(p, ipm.diag, "germ_x", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
               xlab="May Precip", ylab="Germination probability", 
               xlim=range(x2), ylim=c(0,1))
  plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
  legend("center", lty=c(rep(1, length(iss.col)), 3), col=c(iss.COL, 1), cex=1.5,
         lwd=c(rep(3, length(iss.col)+1)), legend=c(names(ipm.diag), "True"))
  if(write.plots) dev.off()
}



########-------------------------------
## CAd: Temperature
########
if(by.issue) {
  for(i in 1:length(issues)) {
    f.name <- paste0("figs/diag/byIssue/", sp_i$Num, "_CAp_Temp_", issues[i], ".jpg")
    if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
    par(mfrow=c(2,3))
    plot_sdm_reg(p, cad.diag[i], "s.M", 1:3, x1, cbind(1, x1.mx), iss.col[i], 
                 xlab="Min Temp: bio6", ylab="Juvenile survival", 
                 xlim=range(x1), ylim=c(0,1))
    plot_sdm_reg(p, cad.diag[i], "K", 1:3, x1, cbind(1, x1.mx), iss.col[i], 
                 xlab="Min Temp: bio6", ylab="Carrying capacity", 
                 xlim=range(x1), ylim=c(0,5e3))
    plot_sdm_reg(p, cad.diag[i], "p.f", 1:3, x1, cbind(1, x1.mx), iss.col[i], 
                 xlab="Min Temp: bio6", ylab="Flowering probability", 
                 xlim=range(x1), ylim=c(0,1))
    plot_sdm_reg(p, cad.diag[i], "mu", 1:3, x1, cbind(1, x1.mx), iss.col[i], 
                 xlab="Min Temp: bio6", ylab="Seed production", 
                 xlim=range(x1), ylim=c(0,1e3))
    plot_sdm_reg(p, cad.diag[i], "g.D", 1:3, x1, cbind(1, x1.mx), iss.col[i], 
                 xlab="Min Temp: bio6", ylab="Germination probability", 
                 xlim=range(x1), ylim=c(0,1))
    plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
    legend("center", lty=rep(1, length(iss.col)), col=iss.COL, cex=1.5,
           lwd=rep(3, length(iss.col)), legend=names(cad.diag))
    if(write.plots) dev.off()
  }
} else {
  f.name <- paste0("figs/diag/", sp_i$Num, "_CAp_Temp.jpg")
  if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
  par(mfrow=c(2,3))
  plot_sdm_reg(p, cad.diag, "s.M", 1:3, x1, cbind(1, x1.mx), iss.col, 
               xlab="Min Temp: bio6", ylab="Juvenile survival", 
               xlim=range(x1), ylim=c(0,1))
  plot_sdm_reg(p, cad.diag, "K", 1:3, x1, cbind(1, x1.mx), iss.col, 
               xlab="Min Temp: bio6", ylab="Carrying capacity", 
               xlim=range(x1), ylim=c(0,5e3))
  plot_sdm_reg(p, cad.diag, "p.f", 1:3, x1, cbind(1, x1.mx), iss.col, 
               xlab="Min Temp: bio6", ylab="Flowering probability", 
               xlim=range(x1), ylim=c(0,1))
  plot_sdm_reg(p, cad.diag, "mu", 1:3, x1, cbind(1, x1.mx), iss.col, 
               xlab="Min Temp: bio6", ylab="Seed production", 
               xlim=range(x1), ylim=c(0,3e2))
  plot_sdm_reg(p, cad.diag, "g.D", 1:3, x1, cbind(1, x1.mx), iss.col, 
               xlab="Min Temp: bio6", ylab="Germination probability", 
               xlim=range(x1), ylim=c(0,1))
  plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
  legend("center", lty=rep(1, length(iss.col)), col=iss.COL, cex=1.5,
         lwd=rep(3, length(iss.col)), legend=names(cad.diag))
  if(write.plots) dev.off()
}



########-------------------------------
## CAd: Precipitation
########
if(by.issue) {
  for(i in 1:length(issues)) {
    f.name <- paste0("figs/diag/byIssue/", sp_i$Num, "_CAp_Precip_", issues[i], ".jpg")
    if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
    par(mfrow=c(2,3))
    plot_sdm_reg(p, cad.diag[i], "s.M", c(1,4,5), x2, cbind(1, x2.mx), iss.col[i], 
                 xlab="May Precip", ylab="Juvenile survival", 
                 xlim=range(x2), ylim=c(0,1))
    plot_sdm_reg(p, cad.diag[i], "K", c(1,4,5), x2, cbind(1, x2.mx), iss.col[i], 
                 xlab="May Precip", ylab="Carrying capacity", 
                 xlim=range(x2), ylim=c(0,5e3))
    plot_sdm_reg(p, cad.diag[i], "p.f", c(1,4,5), x2, cbind(1, x2.mx), iss.col[i], 
                 xlab="May Precip", ylab="Flowering probability", 
                 xlim=range(x2), ylim=c(0,1))
    plot_sdm_reg(p, cad.diag[i], "mu", c(1,4,5), x2, cbind(1, x2.mx), iss.col[i], 
                 xlab="May Precip", ylab="Seed production", 
                 xlim=range(x2), ylim=c(0,6e3))
    plot_sdm_reg(p, cad.diag[i], "g.D", c(1,4,5), x2, cbind(1, x2.mx), iss.col[i], 
                 xlab="May Precip", ylab="Germination probability", 
                 xlim=range(x2), ylim=c(0,1))
    plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
    legend("center", lty=rep(1, length(iss.col)), col=iss.COL, cex=1.5,
           lwd=rep(3, length(iss.col)), legend=names(cad.diag))
    if(write.plots) dev.off()
  }
} else {
  f.name <- paste0("figs/diag/", sp_i$Num, "_CAp_Precip.jpg")
  if(write.plots) jpeg(f.name, width=8, height=6, units="in", res=400)
  par(mfrow=c(2,3))
  plot_sdm_reg(p, cad.diag, "s.M", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
               xlab="May Precip", ylab="Juvenile survival", 
               xlim=range(x2), ylim=c(0,1))
  plot_sdm_reg(p, cad.diag, "K", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
               xlab="May Precip", ylab="Carrying capacity", 
               xlim=range(x2), ylim=c(0,5e3))
  plot_sdm_reg(p, cad.diag, "p.f", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
               xlab="May Precip", ylab="Flowering probability", 
               xlim=range(x2), ylim=c(0,1))
  plot_sdm_reg(p, cad.diag, "mu", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
               xlab="May Precip", ylab="Seed production", 
               xlim=range(x2), ylim=c(0,3e2))
  plot_sdm_reg(p, cad.diag, "g.D", c(1,4,5), x2, cbind(1, x2.mx), iss.col, 
               xlab="May Precip", ylab="Germination probability", 
               xlim=range(x2), ylim=c(0,1))
  plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F, main=sp)
  legend("center", lty=rep(1, length(iss.col)), col=iss.COL, cex=1.5,
         lwd=rep(3, length(iss.col)), legend=names(cad.diag))
  if(write.plots) dev.off()
}





########-------------------------------
## MaxEnt: AUC
########
ggplot(mxe.auc, aes(x=AUC, colour=issue)) + geom_density() +
  facet_wrap(~Boundary) + scale_colour_manual(values=iss.COL) + 
  ggtitle(sp) + xlim(0.95, 1)
ggsave(paste0("figs/diag/", sp_i$Num, "_MxE_AUC.jpg"),
       width=8, height=3, units="in", dpi=400)



########-------------------------------
## MaxEnt: ROC
########
ggplot(mxe.roc, aes(x=FPR, y=TPR, colour=issue, group=obs)) + 
  geom_abline(linetype=3, alpha=0.5) + geom_line(alpha=0.1) +
  facet_grid(Boundary~issue) + scale_colour_manual(values=iss.COL) + 
  labs(main=sp, xlab="1 - specificity", ylab="Sensitivity") + 
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  scale_x_continuous(breaks=c(0,0.5,1)) + scale_y_continuous(breaks=c(0,0.5,1))
ggsave(paste0("figs/diag/", sp_i$Num, "_MxE_ROC.jpg"),
       width=8, height=3, units="in", dpi=400)




