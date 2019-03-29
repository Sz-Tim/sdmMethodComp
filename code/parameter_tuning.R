
sp <- c("barberry", "garlic_mustard")[1]
res <- "5km"
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
# em.df <- st_as_sf(read.csv("data/eddmaps/BETH_ALPE4.csv"), 
#                   coords=c("Longitude_Decimal", "Latitude_Decimal")) %>% 
#   filter(USDAcode==ifelse(sp=="barberry", "BETH", "ALPE4")) 
# em.df <- em.df %>% st_set_crs(4326) %>% 
#   st_transform(., CRS(paste("+proj=aea", "+lat_1=29.5", "+lat_2=45.5", "+lat_0=23",
#                             "+lon_0=-96", "+x_0=0", "+y_0=0", "+ellps=GRS80",
#                             "+towgs84=0,0,0,-0,-0,-0,0", "+units=m", "+no_defs",  
#                             collapse=" "))@projargs) %>%
#   st_coordinates %>% as.data.frame %>% rename(lon=X, lat=Y) %>%
#   mutate(ord=em.df$ord)
# em.id <- lapply(1:nrow(em.df), function(x) get_pt_id(L$env.in, unlist(em.df[x,])))
# em.id <- unlist(em.id)[unlist(em.id)>0]
# em_ids <- unique(em.id)
# em_ord <- map_dbl(em_ids, ~mean(em.df$ord[which(em.id==.)]))
# abund.pred <- dir("data/eddmaps", ".asc", full.names=T) %>% 
  # map(raster::raster) 
hs.df <- read.csv("data/AllenBradley2016/IAS_occurences_final_analysis.csv") %>%
  filter(PLANT_CODE==ifelse(sp=="barberry", "BETH", "ALPE4")) %>% 
  st_as_sf(coords=c("LONGITUDE_DECIMAL", "LATITUDE_DECIMAL")) %>%
  st_set_crs(4326) %>% 
  st_transform(., CRS(paste("+proj=aea", "+lat_1=29.5", "+lat_2=45.5", "+lat_0=23",
                            "+lon_0=-96", "+x_0=0", "+y_0=0", "+ellps=GRS80",
                            "+towgs84=0,0,0,-0,-0,-0,0", "+units=m", "+no_defs",  
                            collapse=" "))@projargs) %>%
  st_coordinates %>% as.data.frame %>% rename(lon=X, lat=Y)
hs.id <- lapply(1:nrow(hs.df), function(x) get_pt_id(L$env.in, unlist(hs.df[x,])))
hs.id <- unlist(hs.id)[unlist(hs.id)>0]
hs_ids <- unique(hs.id)



lam.df <- readRDS(paste0("vs/sp", ifelse(sp=="barberry", 1, 2), "/lam_df.rds"))
p <- readRDS(paste0("vs/sp", ifelse(sp=="barberry", 1, 2), "/p.rds"))
sdd.ji <- readRDS(paste0("vs/sp", ifelse(sp=="barberry", 1, 2), "/sdd_ji.rds"))
p.ji <- readRDS(paste0("vs/sp", ifelse(sp=="barberry", 1, 2), "/p_ji.rds"))


p.pnas <- fit_PNAS_species(sp, env.f, nlcd.sp, clim_X, FALSE, max_z_pow, habitat,
                           x_min, x_max, y_min, y_max)
# p <- p.pnas

if(sp=="garlic_mustard") {
  p$s_x <- c(-2, -0.9, -0.3, -0.6)
  p$g_x <- c(-4, -1.1, -1, -0.4)
  p$germ_x <- c(0.7, -0.5, -0.9, -0.2, -0.1)
  p$fl_x <- c(-0.3, -0.05, 0.25, -0.1)
  p$seed_x <- c(-1.3, -1, -0.3, -0.3)
} else {
  p$s_x <- c(-5, -2.75, 1, -2.5)
  p$g_x <- c(-1.5, -0.4, -0.2, -0.5)
  p$germ_x <- c(-1.25, -4, -2, -2, -0.75)
}

p$n <- 20
n_z <- list(s=length(p$s_z), g=length(p$g_z), 
            fl=length(p$fl_z), seed=length(p$seed_z))
n_x <- list(s=length(p$s_x), g=length(p$g_x), 
            fl=length(p$fl_x), seed=length(p$seed_x), germ=length(p$germ_x))
X <- map(n_x, ~as.matrix(L$env.in[,grep(paste(clim_X, collapse="|"), names(L$env.in))]))
X$germ <- cbind(1, X$germ[,-n_x$germ])

U <- fill_IPM_matrices(n.cell, buffer=0, discrete=1, p, n_z, n_x, 
                       X, sdd.ji, p.ji, verbose=T)
if(sp=="garlic_mustard") {
  library(doSNOW); library(foreach)
  p.c <- makeCluster(4); registerDoSNOW(p.c)
  U$lambda <- foreach(i=1:n.cell, .combine="c") %dopar% {
    iter_lambda(p, U$Ps[,,i], U$Fs[,,i], tol=0.5)
  }
  stopCluster(p.c)
} else {
  U$lambda <- apply(U$IPMs, 3, function(x) Re(eigen(x)$values[1]))
}
#lam.df <- L$env.in %>% mutate(lambda=U$lambda)
lam.df <- lam.df %>% mutate(lambda.new=U$lambda)


# lam.df$em <- FALSE
# lam.df$em[lam.df$id %in% em_ids] <- TRUE
# lam.df$em_ord <- NA
# lam.df$em_ord[match(em_ids, lam.df$id)] <- em_ord
lam.df$hs <- FALSE
lam.df$hs[lam.df$id %in% hs_ids] <- TRUE
lam.df$s <- antilogit(as.matrix(lam.df[,c(11,12,39,40)]) %*% p$s_x)
lam.df$g <- as.matrix(lam.df[,c(11,12,39,40)]) %*% p$g_x
lam.df$germ <- antilogit(cbind(1, as.matrix(lam.df[,c(11,12,39,40)])) %*% p$germ_x)
lam.df$seed <- exp(as.matrix(lam.df[,c(11,12,39,40)]) %*% p$seed_x)
x1 <- seq(min(lam.df$bio10_6), max(lam.df$bio10_6), length.out=200)
x1.mx <- cbind(x1, x1^2)
x2 <- seq(min(lam.df$bio10_prMay), max(lam.df$bio10_prMay), length.out=200)
x2.mx <- cbind(x2, x2^2)

with(filter(lam.df, hs), sum(lambda.new>1)/sum(hs))
ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, lambda.new>1), aes(lon, lat, fill=lambda.new)) +
  scale_fill_viridis(option="B") +
  geom_point(data=filter(lam.df, hs & lambda.new<1), aes(lon, lat), colour="white", shape=1)

summary(filter(lam.df, em)$lambda>1)
ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, lambda>1), aes(lon, lat, fill=lambda)) +
  scale_fill_viridis(option="B") +
  geom_point(data=filter(lam.df, em & lambda<1), aes(lon, lat), colour="white", shape=1)

summary(filter(lam.df, em)$Surv.S>0)
ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, Surv.S>0), aes(lon, lat, fill=Surv.S)) +
  scale_fill_viridis(option="B") +
  geom_point(data=filter(lam.df, em & Surv.S==0), aes(lon, lat), colour="white", shape=1)

summary(filter(lam.df, hs)$Surv.S>0)
ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, Surv.S>0), aes(lon, lat, fill=Surv.S)) +
  scale_fill_viridis(option="B") +
  geom_point(data=filter(lam.df, hs & Surv.S==0), aes(lon, lat), colour="white", shape=1)




par(mfrow=c(2,3))
plot(x1, antilogit(x1.mx %*% p$s_x[1:2]), xlab="Temp", ylab="Survival", type="l",
     ylim=c(0,1))
lines(x1, antilogit(x1.mx %*% p.pnas$s_x[1:2]), col="red")
plot(x1, exp(x1.mx %*% p$seed_x[1:2]), xlab="Temp", ylab="Seeds", type="l")
lines(x1, exp(x1.mx %*% p.pnas$seed_x[1:2]), col="red")
plot(x1, x1.mx %*% p$g_x[1:2], xlab="Temp", ylab="Growth", type="l")
lines(x1, x1.mx %*% p.pnas$g_x[1:2], col="red")
# plot(x1, antilogit(cbind(1, x1.mx) %*% p$germ_x[1:3]), xlab="Temp", 
#      ylab="Germination", type="l", ylim=c(0,1))
# lines(x1, antilogit(cbind(1, x1.mx) %*% p.pnas$germ_x[1:3]), col="red")
plot(x2, antilogit(x2.mx %*% p$s_x[3:4]), xlab="Precip", ylab="Survival", type="l",
     ylim=c(0,1))
lines(x2, antilogit(x2.mx %*% p.pnas$s_x[3:4]), col="red")
plot(x2, exp(x2.mx %*% p$seed_x[3:4]), xlab="Precip", ylab="Seeds", type="l")
lines(x2, exp(x2.mx %*% p.pnas$seed_x[3:4]), col="red")
plot(x2, x2.mx %*% p$g_x[3:4], xlab="Precip", ylab="Growth", type="l")
lines(x2, x2.mx %*% p.pnas$g_x[3:4], col="red")
# plot(x2, antilogit(x2.mx %*% p$germ_x[4:5]), xlab="Precip", 
#      ylab="Germination", type="l", ylim=c(0,1))
# lines(x2, antilogit(x2.mx %*% p.pnas$germ_x[4:5]), col="red")

ggplot(lam.df) + geom_tile(aes(lon, lat, fill=bio10_6)) +
  scale_fill_viridis() +
  geom_point(data=filter(lam.df, hs), 
             aes(lon, lat), colour="white", shape=1)

ggplot(lam.df) + geom_tile(aes(lon, lat, fill=bio10_prMay)) +
  scale_fill_viridis() +
  geom_point(data=filter(lam.df, hs), 
             aes(lon, lat), colour="white", shape=1)

ggplot(lam.df) + geom_tile(aes(lon, lat, fill=s)) +
  scale_fill_viridis(limits=c(0,1)) +
  geom_point(data=filter(lam.df, hs & lambda<1), 
             aes(lon, lat), colour="white", shape=1)
ggplot(lam.df) + geom_tile(aes(lon, lat, fill=g)) +
  scale_fill_viridis() +
  geom_point(data=filter(lam.df, hs & lambda<1), 
             aes(lon, lat), colour="white", shape=1)
ggplot(lam.df) + geom_tile(aes(lon, lat, fill=seed)) +
  scale_fill_viridis() +
  geom_point(data=filter(lam.df, hs & lambda<1), 
             aes(lon, lat), colour="white", shape=1)
ggplot(lam.df) + geom_tile(aes(lon, lat, fill=germ)) +
  scale_fill_viridis(limits=c(0,1)) +
  geom_point(data=filter(lam.df, hs & lambda<1), 
             aes(lon, lat), colour="white", shape=1)


ggplot() + geom_tile(data=lam.df, aes(lon, lat), fill="gray30") +
  geom_tile(data=filter(lam.df, lambda.new>1), aes(lon, lat, fill=lambda.new)) +
  scale_fill_viridis(option="B") + ggtitle("hotspots data") +
  geom_point(data=filter(lam.df, hs), aes(lon, lat), colour="white", 
             shape=1, alpha=0.7)








