
sp <- c("barberry", "garlic_mustard")[1]
res <- "10km"
clim_X <- paste0("bio10_", c(6, "prMay"))
habitat <- 4
max_z_pow <- 1
x_min <- 0#200#675#
x_max <- Inf
y_min <- 0
y_max <- Inf#75#250#

pkgs <- c("gbPopMod", "tidyverse", "rgdal", "here", "sf", "viridis", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "fn", full.names=T), source)
env.f <- paste0("data/ENF_", res, ".csv")
sp_i <- read.csv(paste0("data/species_", res, ".csv")) %>% filter(Name==sp)
nlcd.sp <- read.csv(here("data/PNAS_2017/", sp_i$LC_f))
L <- build_landscape(env.f, nlcd.sp, x_min, x_max, y_min, y_max) 
n.cell <- sum(L$env.rct$inbd)
em.df <- st_as_sf(read.csv("data/eddmaps/BETH_ALPE4.csv"), 
                  coords=c("Longitude_Decimal", "Latitude_Decimal")) %>% 
  filter(USDAcode==ifelse(sp=="barberry", "BETH", "ALPE4")) 
em.df <- em.df %>% st_set_crs(4326) %>% 
  st_transform(., CRS(paste("+proj=aea", "+lat_1=29.5", "+lat_2=45.5", "+lat_0=23",
                            "+lon_0=-96", "+x_0=0", "+y_0=0", "+ellps=GRS80",
                            "+towgs84=0,0,0,-0,-0,-0,0", "+units=m", "+no_defs",  
                            collapse=" "))@projargs) %>%
  st_coordinates %>% as.data.frame %>% rename(lon=X, lat=Y) %>%
  mutate(ord=em.df$ord)
  
em.id <- lapply(1:nrow(em.df), function(x) get_pt_id(L$env.in, unlist(em.df[x,])))
em.id <- unlist(em.id)[unlist(em.id)>0]
em_ids <- unique(em.id)
em_ord <- map_dbl(em_ids, ~mean(em.df$ord[which(em.id==.)]))
abund.pred <- dir("data/eddmaps", ".asc", full.names=T) %>% 
  map(raster::raster) 
  


p.pnas <- fit_PNAS_species(sp, env.f, nlcd.sp, clim_X, FALSE, max_z_pow, habitat,
                           x_min, x_max, y_min, y_max)
p <- p.pnas
# too high in the south
p$s_x <- c(p$s_x[1], p$s_x[2], p$s_x[3], p$s_x[4])
p$g_x <- c(p$g_x[1], p$g_x[2], p$g_x[3], p$g_x[4])
p$germ_x <- c(p$germ_x[1], p$germ_x[2], p$germ_x[3], p$germ_x[4], p$germ_x[5])

p$s_x <- c(-2.75, -1.25, 3.5, 0)
p$g_x <- c(-1.25, -0.75, -0.5, -0.3)
p$germ_x <- c(1.5, -4, -1.75, -2.5, 0)

p$n <- 20
p$p_emig <- 0
n_z <- list(s=length(p$s_z), g=length(p$g_z), 
            fl=length(p$fl_z), seed=length(p$seed_z))
n_x <- list(s=length(p$s_x), g=length(p$g_x), 
            fl=length(p$fl_x), seed=length(p$seed_x), germ=length(p$germ_x))
X <- map(n_x, ~as.matrix(L$env.in[,grep(paste(clim_X, collapse="|"), names(L$env.in))]))
if(!is.null(X$germ)) X$germ <- cbind(1, X$germ[,-n_x$germ])

U <- fill_IPM_matrices(n.cell, buffer=0, discrete=1, p, n_z, n_x, 
                       X, sdd.ji, p.ji, sp, verbose=T)
if(sp=="garlic_mustard") {
  U$lambda <- mclapply(1:n.cell, 
                       function(x) iter_lambda(p, U$Ps[,,x], U$Fs[,,x]),
                       mc.cores=4) %>% simplify2array()
} else {
  U$lambda <- apply(U$IPMs, 3, function(x) Re(eigen(x)$values[1]))
}
lam.df <- L$env.in %>% mutate(lambda=U$lambda)
lam.df$em <- FALSE
lam.df$em[lam.df$id %in% em_ids] <- TRUE
lam.df$em_ord <- NA
lam.df$em_ord[match(em_ids, lam.df$id)] <- em_ord
x1 <- seq(min(lam.df$bio10_6), max(lam.df$bio10_6), length.out=200)
x1.mx <- cbind(x1, x1^2)
x2 <- seq(min(lam.df$bio10_prMay), max(lam.df$bio10_prMay), length.out=200)
x2.mx <- cbind(x2, x2^2)

summary(filter(lam.df, em)$lambda>1)
ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, lambda>1), aes(lon, lat, fill=lambda)) +
  scale_fill_viridis(option="B") +
  geom_point(data=filter(lam.df, em & lambda<1), aes(lon, lat), colour="white", shape=1)
par(mfrow=c(1,1)); plot(abund.pred[[3]], main="Barberry")

par(mfrow=c(2,3))
plot(x1, antilogit(x1.mx %*% p$s_x[1:2]), xlab="Temp", ylab="Survival", type="l",
     ylim=c(0,1))
lines(x1, antilogit(x1.mx %*% p.pnas$s_x[1:2]), col="red")
plot(x1, x1.mx %*% p$g_x[1:2], xlab="Temp", ylab="Growth", type="l")
lines(x1, x1.mx %*% p.pnas$g_x[1:2], col="red")
plot(x1, antilogit(cbind(1, x1.mx) %*% p$germ_x[1:3]), xlab="Temp", 
     ylab="Germination", type="l", ylim=c(0,1))
lines(x1, antilogit(cbind(1, x1.mx) %*% p.pnas$germ_x[1:3]), col="red")
plot(x2, antilogit(x2.mx %*% p$s_x[3:4]), xlab="Precip", ylab="Survival", type="l",
     ylim=c(0,1))
lines(x2, antilogit(x2.mx %*% p.pnas$s_x[3:4]), col="red")
plot(x2, x2.mx %*% p$g_x[3:4], xlab="Precip", ylab="Growth", type="l")
lines(x2, x2.mx %*% p.pnas$g_x[3:4], col="red")
plot(x2, antilogit(x2.mx %*% p$germ_x[4:5]), xlab="Precip", 
     ylab="Germination", type="l", ylim=c(0,1))
lines(x2, antilogit(x2.mx %*% p.pnas$germ_x[4:5]), col="red")

ggplot(lam.df) + geom_tile(aes(lon, lat, fill=bio10_6)) +
  scale_fill_viridis()
ggplot(lam.df) + geom_tile(aes(lon, lat, fill=bio10_prMay)) +
  scale_fill_viridis()



