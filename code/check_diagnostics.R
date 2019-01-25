
sp <- c("barberry", "garlic_mustard")[1]
pkgs <- c("tidyverse", "magrittr", "stringr", "here", "viridis")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
sp_i <- read.csv("data/species_3km.csv") %>% filter(Name==sp)
vs.dir <- paste0("vs/", sp_i$Num, "/")
p <- readRDS(here(vs.dir, "p.rds"))
lam.df <- readRDS(here(vs.dir, "lam_df.rds"))

ipm.diag <- setNames(dir("out/sp1", "Diag_IPM", full.names=T) %>% map(readRDS),
                     str_remove(str_remove(dir("out/sp1", "Diag_IPM"), ".rds"), "Diag_IPM_"))
ipm.cols <- RColorBrewer::brewer.pal(length(ipm.diag), "Dark2")
ipm.cols <- apply(col2rgb(ipm.cols)/255, 2, function(x) rgb(x[1],x[2],x[3],0.5)) 

z.seq <- seq(p$z.rng[1], p$z.rng[2], length.out=200)
z.mx <- cbind(1, z.seq)
x <- seq(min(lam.df$bio10_5), max(lam.df$bio10_5), length.out=200)
x.mx <- cbind(x, x^2)

########-------------------------------
## Size
########
par(mfrow=c(2,3))
# survival vs size
plot(NA, NA, xlab="Size (t)", ylab="Survival", ylim=c(0,1), xlim=p$z.rng, cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(z.seq, antilogit(z.mx %*% .$s_z), col=ipm.cols[i]))
}
lines(z.seq, antilogit(z.mx %*% p$s_z), lwd=2, lty=3)
  
# growth vs size
plot(NA, NA, xlab="Size (t)", ylab="Size (t+1)", ylim=p$z.rng, xlim=p$z.rng, cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(z.seq, z.mx %*% .$g_z, col=ipm.cols[i]))
}
lines(z.seq, z.mx %*% p$g_z, lwd=2, lty=3)

# flowering vs size
plot(NA, NA, xlab="Size (t)", ylab="Flowering", ylim=c(0,1), xlim=p$z.rng, cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(z.seq, antilogit(z.mx %*% .$fl_z), col=ipm.cols[i]), cex.lab=1.5)
}
lines(z.seq, antilogit(z.mx %*% p$fl_z), lwd=2, lty=3)

# seeds vs size
plot(NA, NA, xlab="Size (t)", ylab="Seed production", ylim=c(0,6e3), xlim=p$z.rng, cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(z.seq, exp(z.mx %*% .$seed_z), col=ipm.cols[i]))
}
lines(z.seq, exp(z.mx %*% p$seed_z), lwd=2, lty=3)

# recruit size distribution
plot(NA, NA, xlab="Recruit size", ylab="density", xlim=p$z.rng, ylim=c(0, 0.4), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~curve(dnorm(x, .$rcr_z[1], .$rcr_z[2]), add=T, 
                             from=p$z.rng[1], to=p$z.rng[2], col=ipm.cols[i]))
}
curve(dnorm(x, p$rcr_z[1], p$rcr_z[2]), from=p$z.rng[1], to=p$z.rng[2], lwd=2, 
      add=T, lty=3)

# legend
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=c(rep(1, length(ipm.cols)), 3), lwd=c(rep(3, length(ipm.cols)+1)),
       col=c(ipm.cols, 1), legend=c(names(ipm.diag), "True"), cex=1.5)



########-------------------------------
## Temperature
########
par(mfrow=c(2,3))
# survival vs temp
plot(NA, NA, xlab="Temperature", ylab="Survival", ylim=c(0,1), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, antilogit(x.mx %*% .$s_x[1:2]), col=ipm.cols[i]))
}
lines(x, antilogit(x.mx %*% p$s_x[1:2]), lwd=2, lty=3)

# growth vs temp
plot(NA, NA, xlab="Temperature", ylab="Growth", ylim=c(-2,2), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, x.mx %*% .$g_x[1:2], col=ipm.cols[i]))
}
lines(x, x.mx %*% p$g_x[1:2], lwd=2, lty=3)

# flowering vs temp
plot(NA, NA, xlab="Temperature", ylab="Flowering", ylim=c(0,1), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, antilogit(x.mx %*% .$fl_x[1:2]), col=ipm.cols[i]))
}
lines(x, antilogit(x.mx %*% p$fl_x[1:2]), lwd=2, lty=3)

# seeds vs temp
plot(NA, NA, xlab="Temperature", ylab="Seed production", ylim=c(-2,2), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, exp(x.mx %*% .$seed_x[1:2]), col=ipm.cols[i]))
}
lines(x, exp(x.mx %*% p$seed_x[1:2]), lwd=2, lty=3)

# germination vs temp
plot(NA, NA, xlab="Temperature", ylab="Germination", ylim=c(0,1), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, antilogit(cbind(1, x.mx) %*% .$germ_x[1:3]), col=ipm.cols[i]))
}
lines(x, antilogit(cbind(1, x.mx) %*% p$germ_x[1:3]), lwd=2, lty=3)

# legend
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=c(rep(1, length(ipm.cols)), 3), lwd=c(rep(3, length(ipm.cols)+1)),
       col=c(ipm.cols, 1), legend=c(names(ipm.diag), "True"), cex=1.5)





########-------------------------------
## Precipitation
########
par(mfrow=c(2,3))
# survival vs temp
plot(NA, NA, xlab="Precipitation", ylab="Survival", ylim=c(0,1), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, antilogit(x.mx %*% .$s_x[3:4]), col=ipm.cols[i]))
}
lines(x, antilogit(x.mx %*% p$s_x[3:4]), lwd=2, lty=3)

# growth vs temp
plot(NA, NA, xlab="Precipitation", ylab="Growth", ylim=c(-10,5), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, x.mx %*% .$g_x[3:4], col=ipm.cols[i]))
}
lines(x, x.mx %*% p$g_x[3:4], lwd=2, lty=3)

# flowering vs temp
plot(NA, NA, xlab="Precipitation", ylab="Flowering", ylim=c(0,1), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, antilogit(x.mx %*% .$fl_x[3:4]), col=ipm.cols[i]))
}
lines(x, antilogit(x.mx %*% p$fl_x[3:4]), lwd=2, lty=3)

# seeds vs temp
plot(NA, NA, xlab="Precipitation", ylab="Seed production", ylim=c(-2,2), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, exp(x.mx %*% .$seed_x[3:4]), col=ipm.cols[i]))
}
lines(x, exp(x.mx %*% p$seed_x[3:4]), lwd=2, lty=3)

# germination vs temp
plot(NA, NA, xlab="Precipitation", ylab="Germination", ylim=c(0,1), xlim=range(x), cex.lab=1.5)
for(i in seq_along(ipm.diag)) {
  walk(ipm.diag[[i]], ~lines(x, antilogit(cbind(1, x.mx) %*% .$germ_x[c(1,4:5)]), col=ipm.cols[i]))
}
lines(x, antilogit(cbind(1, x.mx) %*% p$germ_x[c(1,4:5)]), lwd=2, lty=3)

# legend
plot(NA, NA, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1), axes=F)
legend("center", lty=c(rep(1, length(ipm.cols)), 3), lwd=c(rep(3, length(ipm.cols)+1)),
       col=c(ipm.cols, 1), legend=c(names(ipm.diag), "True"), cex=1.5)

