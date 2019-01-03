# NSF IIS 
# GIS processing
# Tim Szewczyk


library(maptools); library(raster); library(rgeos); library(rgdal)
library(spatialEco); library(tidyverse); library(doSNOW); library(foreach)
n_core <- 8
cell_side <- 3000  # cell dimensions in meters; cell_area = cell_side^2



##------ Map details
## Assign CRS definitions, load and project state map
##------
# Albers equal area definition from NLCD
alb_CRS <- CRS(paste("+proj=aea", "+lat_1=29.5", "+lat_2=45.5", "+lat_0=23",
                 "+lon_0=-96", "+x_0=0", "+y_0=0", "+ellps=GRS80",
                 "+towgs84=0,0,0,-0,-0,-0,0", "+units=m", "+no_defs",  
                 collapse=" "))
# WGS84 definition
wgs_CRS <- CRS("+init=epsg:4326")
# Map of state boundaries projected to alb_CRS
states <- maps::map("state", plot=F, fill=T)
IDs <- states$names
states <- map2SpatialPolygons(states, proj4string=wgs_CRS, IDs=IDs)
states <- spTransform(states, alb_CRS)
states <- gBuffer(states, byid=T, width=0)



##------ Ecoregion: Define extent
## Data for ecoregions contained in data/ecoregions/NA_CEC_Eco_Level1.shp. Files
## were downloaded on 2018 June 11 from:
## <ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l1.zip>
##------
# All ecoregions in North America
Ecoregions <- readOGR(dsn="data/ecoregions", layer="NA_CEC_Eco_Level1") %>%
  spTransform(alb_CRS)
# Eastern Temperate Forest & Northern Forest
ENF <- subset(Ecoregions, NA_L1CODE %in% c(5,8))
ENF.us <- gIntersection(states, ENF, byid=F)
ENF.wgs <- spTransform(ENF.us, wgs_CRS)
ENF.ext <- extent(as.vector(t(bbox(ENF.us))))
ENF.box <- as(ENF.ext, "SpatialPolygons")



##------ Initialize grid
## Grid resolution in meters, extent is ENF rectangle
##------
grd <- raster(ext=ENF.ext, crs=alb_CRS, resolution=cell_side) %>%
  rasterize(ENF.us, .) %>% 
  rasterToPolygons(., n=4)
grd@data$layer <- 1:ncell(grd)
grd.wgs <- spTransform(grd, wgs_CRS)
grd.df <- as.data.frame(coordinates(grd))
names(grd.df) <- c("long", "lat")

clim.rast <- raster("data/climate/CHELSA_prec_5_V1.2_land.tif")
# run crop(), projectRaster() below
clim.df <- zonal.stats(grd, clim.rast, mean, trace=F, plot=F)

##------ Climate data
## Climate raster files are stored in climate/ with no other files. Each file
## is named CHELSA_bio10_X.tif. If names are different, adjust arguments in
## str_split_fixed(). Bioclimatic variables were downloaded on 2018 June 11 
## from CHELSA <http://chelsa-climate.org/downloads/>
##------
# Load and extract bioclimatic variable numbers
clim.f <- dir("data/climate", "bio10")
clim.var <- str_split_fixed(clim.f, pattern="[[:punct:]]", 4)[,2:3] %>% 
  apply(1, str_flatten, collapse="_")
clim.rast <- stack(paste0("data/climate/", clim.f))
# Crop and project
clim.rast <- crop(clim.rast, extent(ENF.wgs), snap="out")
clim.rast <- projectRaster(clim.rast, crs=alb_CRS)
# Calculate mean value within each grid cell
p.c <- makeCluster(n_core); registerDoSNOW(p.c)
clim.df <- foreach(i=seq_along(clim.var), .combine="data.frame", .errorhandling="pass",
                   .packages=c("sp", "spatialEco")) %dopar% {
  zonal.stats(grd, clim.rast[[i]], mean, trace=F, plot=F)
}
stopCluster(p.c)
names(clim.df) <- clim.var
grd.df <- cbind(grd.df, clim.df)
rm(clim.rast); removeTmpFiles(0)



##------ Landcover data
## NLCD files are stored in data/nlcd/. The NLCD 2011 Land Cover files were 
## downloaded on 2018 June 11 from <https://www.mrlc.gov/nlcd11_data.php>. If
## aggregating land cover categories first, the script requires a key in
## data/nlcd/reclass.txt where the first column is the NLCD 2011 numeric code
## and the second column is the aggregated category number.
##------
# Load; already projected in alb_CRS
nlcd.rast <- raster("data/nlcd/nlcd_2011_landcover_2011_edition_2014_10_10.img")
nlcd.rast <- crop(nlcd.rast, extent(ENF.us), snap="out")
# For reclassifying before summarizing
reclass_first <- FALSE
if(reclass_first) {
  nlcd.reclass <- as.matrix(read.table("data/nlcd/reclass.txt", header=T))
  nlcd.rast <- reclassify(nlcd.rast, nlcd.reclass)
  nlcd.val <- sort(unique(nlcd.reclass[,2]))
} else {
  # Extract unique values
  nlcd.val <- dplyr::filter(nlcd.rast@data@attributes[[1]], COUNT>0 & ID!=0)$ID
}
# Separate land cover categories to binary layers & calculate mean (=proportion)
p.c <- makeCluster(n_core); registerDoSNOW(p.c)
nlcd.df <- foreach(i=seq_along(nlcd.val), .errorhandling="pass",
                   .packages=c("raster", "sp", "spatialEco")) %dopar% {
  f.i <- paste0("data/nlcd/layers/lc_", nlcd.val[i], ".grd")
  if(file.exists(f.i)) {
    nlcd.i <- raster(f.i)
  } else {
    nlcd.i <- layerize(nlcd.rast, classes=nlcd.val[i], filename=f.i)
  }
  zonal.stats(grd, nlcd.i, mean, trace=F, plot=F)
}
stopCluster(p.c)
rm(nlcd.rast)
nlcd.df <- as.data.frame(do.call("cbind", nlcd.df)) 
nlcd.df <- nlcd.df/rowSums(nlcd.df)
names(nlcd.df) <- paste0("nlcd_", nlcd.val)
grd.df <- cbind(grd.df, nlcd.df)
removeTmpFiles(0)



##------ Road data
## 2017 Tiger/Line Road shapefiles are stored in data/roads/. The files contain
## lines for all primary and secondary roads by state and were downloaded on 
## 2018 June 13 from <ftp://ftp2.census.gov/geo/tiger/TIGER2017/PRISECROADS>.
##------
library(sf)
road.f <- dir("data/roads", "*.shp$", full.names=T)
rd.grd <- lapply(road.f, st_read) %>%
  do.call(rbind, .) %>%
  st_transform(., alb_CRS@projargs) %>%
  st_crop(., extent(ENF.us)) %>%
  st_intersection(., st_as_sf(grd)) %>%
  mutate(length=st_length(.)) %>%
  group_by(layer) %>%
  summarise(length=sum(length))
grd.df$rd_len <- 0
grd.df$rd_len[rd.grd$layer] <- rd.grd$length



write_csv(grd.df, paste0("data/ENF_", cell_side/1000, "km_pr_may.csv"))





