# PNAS data

########
## Setup
########
# Load packages, plot information, and data. 

pkgs <- c("tidyverse", "magrittr", "here", "rgdal", "sf")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(paste0("code/fn_", c("IPM", "aux", "sim"), ".R"), ~source(here(.)))
sp.f <- dir("data/PNAS_2017", "data_frame", full.names=T)
sp.ls <- map(sp.f, read_csv)
names(sp.ls) <- str_split_fixed(dir("data/PNAS_2017", "data_frame"), 
                                "_data_frame_", 2)[,1]
map(sp.ls, ~sort(unique(.$Plot)))
map(sp.ls, ~sort(unique(.$Location.x)))
map(sp.ls, dim)




########
## Clean data
########
# Location.x (i.e., regions) are not always correct or standardized, and 
# some plots are either NA or do not match one of the 21 plots listed in the
# appendix, and so there are no coordinates. This chunk corrects the regions
# based on Table A.4 in the appendix, and removes rows where the plot is either
# NA or an unknown plot.

plots <- list(WCT=c("HILL", "PGA", "PGB", "RED", "RT7A", "RT7B"),
              CT=c("BMC", "BMD", "JC", "JD", "PSA", "PSB", "SMFC"),
              ME=c("AFA", "AFB", "LYA", "LYB"),
              VT=c("HORTC", "HORTD", "JERC", "JERD"))
for(i in seq_along(sp.ls)) {
  # correct plot names
  for(j in seq_along(plots)) {
    sp.ls[[i]]$Location.x[sp.ls[[i]]$Plot %in% plots[[j]]] <- names(plots)[j]
  }
  # remove rows with unknown plots
  sp.ls[[i]] <- filter(sp.ls[[i]], Plot %in% unlist(plots))
}
map(sp.ls, ~sort(unique(.$Plot)))
map(sp.ls, ~sort(unique(.$Location.x)))
map(sp.ls, dim)
saveRDS(sp.ls, "data/PNAS_2017/species_data_list.rds")




########
## Plot coordinates
########
# Determine lat/long coordinates for each plot. Plot info is stored as an sf 
# object projected to an Albers Equal Area projection to match the projection
# of NLCD.
alb_CRS <- CRS(paste("+proj=aea", "+lat_1=29.5", "+lat_2=45.5", "+lat_0=23",
                     "+lon_0=-96", "+x_0=0", "+y_0=0", "+ellps=GRS80",
                     "+towgs84=0,0,0,-0,-0,-0,0", "+units=m", "+no_defs",  
                     collapse=" ")) # Albers Equal Area from NLCD
plot_i <- sp.ls[[3]] %>% 
  group_by(Location.x, Plot, wplot, lat, long, habitat) %>% 
  summarise %>% ungroup %>%
  mutate(habitat=factor(habitat, labels=c("closed", "open"))) %>%
  rename(canopy=habitat, region=Location.x) %>%
  st_as_sf(., coords=c("long", "lat"), crs=4326) %>%
  st_transform(., crs=alb_CRS@projargs) %>%
  cbind(., st_coordinates(.)) %>%
  rename(lon=X, lat=Y)
write_csv(st_set_geometry(plot_i, NULL), "data/PNAS_2017/plot_coords.csv")



