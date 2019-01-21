ipm.diag <- readRDS("out/sp1/Diag_IPM_none.rds")

z.seq <- seq(p$z.rng[1], p$z.rng[2], length.out=200)
z.mx <- cbind(1, z.seq)
x <- seq(min(lam.df$bio10_5), max(lam.df$bio10_5), length.out=200)
x.mx <- cbind(x, x^2)

########-------------------------------
## Size
########
par(mfrow=c(2,3))
# survival vs size
plot(NA, NA, xlab="Size", ylab="Survival", ylim=c(0,1), xlim=p$z.rng)
ipm.diag %>% walk(~lines(z.seq, antilogit(z.mx %*% .[[1]]$s_z), col="red"))
lines(z.seq, antilogit(z.mx %*% p$s_z), lwd=2)
  
# growth vs size
plot(NA, NA, xlab="Size", ylab="Size Next", ylim=p$z.rng, xlim=p$z.rng)
ipm.diag %>% walk(~lines(z.seq, z.mx %*% .[[1]]$g_z, col="red"))
lines(z.seq, z.mx %*% p$g_z, lwd=2)

# flowering vs size
plot(NA, NA, xlab="Size", ylab="Flowering", ylim=c(0,1), xlim=p$z.rng)
ipm.diag %>% walk(~lines(z.seq, antilogit(z.mx %*% .[[1]]$fl_z), col="red"))
lines(z.seq, antilogit(z.mx %*% p$fl_z), lwd=2)

# seeds vs size
plot(NA, NA, xlab="Size", ylab="Seed production", ylim=c(0,6e3), xlim=p$z.rng)
ipm.diag %>% walk(~lines(z.seq, exp(z.mx %*% .[[1]]$seed_z), col="red"))
lines(z.seq, exp(z.mx %*% p$seed_z), lwd=2)

# recruit size distribution
plot(NA, NA, xlab="Recruit size", ylab="density", xlim=p$z.rng, ylim=c(0, 0.5))
ipm.diag %>% walk(~curve(dnorm(x, .[[1]]$rcr_z[1], .[[1]]$rcr_z[2]), 
                         from=p$z.rng[1], to=p$z.rng[2], add=T, col="red"))
curve(dnorm(x, p$rcr_z[1], p$rcr_z[2]), from=p$z.rng[1], to=p$z.rng[2], lwd=2, add=T)



########-------------------------------
## Temperature
########
par(mfrow=c(2,3))
# survival vs temp
plot(NA, NA, xlab="Temperature", ylab="Survival", ylim=c(0,1), xlim=range(x))
ipm.diag %>% walk(~lines(x, antilogit(x.mx %*% .[[1]]$s_x[1:2]), col="red"))
lines(x, antilogit(x.mx %*% p$s_x[1:2]), lwd=2)

# growth vs temp
plot(NA, NA, xlab="Temperature", ylab="Growth", ylim=c(-2,2), xlim=range(x))
ipm.diag %>% walk(~lines(x, x.mx %*% .[[1]]$g_x[1:2], col="red"))
lines(x, x.mx %*% p$g_x[1:2], lwd=2)

# flowering vs temp
plot(NA, NA, xlab="Temperature", ylab="Flowering", ylim=c(0,1), xlim=range(x))
ipm.diag %>% walk(~lines(x, antilogit(x.mx %*% .[[1]]$fl_x[1:2]), col="red"))
lines(x, antilogit(x.mx %*% p$fl_x[1:2]), lwd=2)

# seeds vs temp
plot(NA, NA, xlab="Temperature", ylab="Seed production", ylim=c(-2,2), xlim=range(x))
ipm.diag %>% walk(~lines(x, exp(x.mx %*% .[[1]]$seed_x[1:2]), col="red"))
lines(x, exp(x.mx %*% p$seed_x[1:2]), lwd=2)

# germination vs temp
plot(NA, NA, xlab="Temperature", ylab="Germination", ylim=c(0,1), xlim=range(x))
ipm.diag %>% walk(~lines(x, antilogit(cbind(1, x.mx) %*% .[[1]]$germ_x[1:3]), col="red"))
lines(x, antilogit(cbind(1, x.mx) %*% p$germ_x[1:3]), lwd=2)




########-------------------------------
## Precipitation
########
par(mfrow=c(2,3))
# survival vs temp
plot(NA, NA, xlab="Precipitation", ylab="Survival", ylim=c(0,1), xlim=range(x))
ipm.diag %>% walk(~lines(x, antilogit(x.mx %*% .[[1]]$s_x[3:4]), col="red"))
lines(x, antilogit(x.mx %*% p$s_x[3:4]), lwd=2)

# growth vs temp
plot(NA, NA, xlab="Precipitation", ylab="Growth", ylim=c(-5,5), xlim=range(x))
ipm.diag %>% walk(~lines(x, x.mx %*% .[[1]]$g_x[3:4], col="red"))
lines(x, x.mx %*% p$g_x[3:4], lwd=2)

# flowering vs temp
plot(NA, NA, xlab="Precipitation", ylab="Flowering", ylim=c(0,1), xlim=range(x))
ipm.diag %>% walk(~lines(x, antilogit(x.mx %*% .[[1]]$fl_x[3:4]), col="red"))
lines(x, antilogit(x.mx %*% p$fl_x[3:4]), lwd=2)

# seeds vs temp
plot(NA, NA, xlab="Precipitation", ylab="Seed production", ylim=c(-2,2), xlim=range(x))
ipm.diag %>% walk(~lines(x, exp(x.mx %*% .[[1]]$seed_x[3:4]), col="red"))
lines(x, exp(x.mx %*% p$seed_x[3:4]), lwd=2)

# germination vs temp
plot(NA, NA, xlab="Precipitation", ylab="Germination", ylim=c(0,1), xlim=range(x))
ipm.diag %>% walk(~lines(x, antilogit(cbind(1, x.mx) %*% .[[1]]$germ_x[c(1,4:5)]), col="red"))
lines(x, antilogit(cbind(1, x.mx) %*% p$germ_x[c(1,4:5)]), lwd=2)

