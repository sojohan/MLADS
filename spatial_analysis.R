#install.packages("spatstat")
library(spatstat)
#install.packages("ggmap")
library(ggmap)
#install.packages("GISTools")
library(GISTools)
#install.packages("sp")
library(sp)
#Til Spatstat

library(rgdal)
library(spdep)
library(gstat)
library(automap)
require(dplyr)


# Spatial Kernel Density Estimation
violent_crimes <- subset(crime, offense != "auto theft" & offense != "theft" & offense != "burglary") 

# order violent crimes 
violent_crimes$offense <- factor(violent_crimes$offense, levels = c("robbery", "aggravated assault", "rape", "murder"))
# restrict to downtown  
violent_crimes <- subset(violent_crimes, -95.39681 <= lon & lon <= -95.34188 & 29.73631 <= lat & lat <= 29.78400)

theme_set(theme_bw(16)) 
HoustonMap <- qmap("houston", zoom = 14, color = "bw", legend = "topleft",source="osm")
HoustonMap <- qmap("houston", zoom = 14, color = "bw", source = "osm", legend="topleft")

get_openstreetmap("houston")
# Density maps From presentation
HoustonMap + geom_point(aes(x = lon, y = lat, colour = offense, size = offense), data = violent_crimes)
HoustonMap + stat_bin2d( aes(x = lon, y = lat, colour = offense, fill = offense), size = .5, bins = 30, alpha = 1/2, data = violent_crimes )
houston <- get_map("houston", zoom = 14) 
HoustonMap <- ggmap("houston", extent = "device", legend = "topleft")
HoustonMap + stat_density2d( aes(x = lon, y = lat, fill = ..level.., alpha = ..level..), size = 2, bins = 4, data = violent_crimes, geom = "polygon" )


#over days
houston <- get_map(location = "houston", zoom = 14, color = "bw", source = "osm")
HoustonMap <- ggmap(houston, base_layer = ggplot(aes(x = lon, y = lat), data = violent_crimes))
HoustonMap + stat_density2d(aes(x = lon, y = lat, fill = ..level.., alpha = ..level..), bins = 5, geom = "polygon", data = violent_crimes) + scale_fill_gradient(low = "black", high = "red") + facet_wrap(~ day)

#data(crime)

maxlon<-max(violent_crimes$lon,na.rm=TRUE)
minlon<-min(violent_crimes$lon, na.rm=TRUE)
maxlat<-max(violent_crimes$lat, na.rm=TRUE)
minlat<-min(violent_crimes$lat, na.rm = TRUE)
na_loc<-complete.cases(violent_crimes$lon)
clean_crime<-violent_crimes[na_loc==TRUE,]

clean_crime<-select(clean_crime,lon,lat,offense)

crimedata<-as.ppp(clean_crime, c(minlon,maxlon, minlat, maxlat ))
crimedata <- crimedata[!duplicated(crimedata)]
plot(crimedata, chars = crimedata$offense, cols = c("red", "blue", "green", "black"))
summary(crimedata)
# Distance to Nearest Event Section
# Density use of spatstat package
plot(density(crimedata))
plot(relrisk(crimedata))
# The K function spatstat package
plot(Kest(crimedata, correction = "border"))
#kcrime<-Kest(crimedata, correction = "border")
kf_crime<-envelope(crimedata,Kest, correction="border")
plot(kf_crime)
mad.test(crimedata,Kest)

# California housing data
# This the data from "The Elements of Statistical learning by Tibshirani, Hastie, Friedman" page 371
 
houses <- read.table("c:/projects/mlads/houses_data.txt", header = FALSE, sep = "",
   dec = ".", row.names = NULL,
   col.names = c("value", "income", "age", "rooms", "bedrooms",
   "pop", "hh", "latitude", "longitude"))


# Interpolation part
# Plot of Californa house prize data
f<-ggplot(data=houses, aes(longitude,latitude,colour=value))
f+geom_point(aes(longitude,latitude))+coord_fixed(1.3)

# Observation box for ppp
maxlon1 <- max(houses$longitude, na.rm = TRUE)
minlon1 <- min(houses$longitude, na.rm = TRUE)
maxlat1 <- max(houses$latitude, na.rm = TRUE)
minlat1 <- min(houses$latitude, na.rm = TRUE)
clean_california <- select(houses, longitude, latitude, value)
caldata <- as.ppp(clean_california, c(minlon1, maxlon1, minlat1, maxlat1))
caldata <- caldata[!duplicated(caldata)]
summary(caldata)
# Idw interpolation calculation
plot(spatstat::idw(caldata))


# Meuse data for kriging. The data is from the sp package
# Krige Part
data(meuse)
coordinates(meuse) = ~x + y
data(meuse.grid)
gridded(meuse.grid) = ~x + y
m <- vgm(.59, "Sph", 874, .04)
# ordinary kriging:
x <- krige(log(zinc) ~ 1, meuse, meuse.grid, model = m)
spplot(x["var1.pred"], main = "ordinary kriging predictions")
# idw part
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE
zn.idw <- krige(log(zinc) ~ 1, meuse, meuse.grid)
library(RColorBrewer)
data(meuse.riv)
meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
meuse.pol <- SpatialPolygons(meuse.lst)
grays <- gray.colors(4, 0.55, 0.95)
cols <- brewer.pal(4, "Reds")
image(zn.idw, col = cols, breaks = log(c(100, 200, 400, 800, 1800)))
plot(meuse.pol, add = TRUE)
plot(meuse, pch = 1, cex = sqrt(meuse$zinc) / 20, add = TRUE)
legVals <- c(100, 200, 500, 1000, 2000)
legend("left", legend = legVals, pch = 1, pt.cex = sqrt(legVals) / 20, bty = "n", title = "measured, ppm", cex = 0.8, y.inter = 0.5)
legend("topleft", fill = cols, legend = c("100-200", "200-400", "400-800", "800-1800"), bty = "n", title = "interpolated, ppm", cex = 0.8, y.inter = 0.5)
title("measured and interpolated zinc")

# Spatial Prediction Models section
setwd("c:/projects/mlads/lat_bundle")
NY8 <- readOGR(".", "NY8_utm18")
TCE <- readOGR(".", "TCE")
NY_nb <- read.gal("NY_nb.gal", region.id = row.names(NY8))
cities <- readOGR(".", "NY8cities")

plot(NY8, border = "grey60", axes = TRUE)
text(coordinates(cities), labels = as.character(cities$names), font = 2, cex = 0.9)
text(bbox(NY8)[1, 1], bbox(NY8)[2, 2], labels = "a)", cex = 0.8)
plot(NY8, border = "grey60", axes = TRUE)
points(TCE, pch = 1, cex = 0.7)
points(TCE, pch = 3, cex = 0.7)
text(coordinates(TCE), labels = as.character(TCE$name), cex = 0.7,
 font = 1, pos = c(4, 1, 4, 1, 4, 4, 4, 2, 3, 4, 2), offset = 0.3)
text(bbox(NY8)[1, 1], bbox(NY8)[2, 2], labels = "b)", cex = 0.8)


plot(NY8, border = "grey60", axes = TRUE)
plot(NY_nb, coordinates(NY8), pch = 19, cex = 0.6, add = TRUE)

# moran test
moran.test(NY8$Cases, listw = nb2listw(NY_nb))

# SAR model
NYlistw <- nb2listw(NY_nb, style = "B")
nysar <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = NY8, listw = NYlistw)
summary(nysar)

library(RColorBrewer)

NY8$sar_trend <- nysar$fit$signal_trend
NY8$sar_stochastic <- nysar$fit$signal_stochastic
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
tr_at <- seq(-1, 1.3, length.out = 21)
tr_rds <- rds(sum(tr_at >= 0) * 2)[ - (1:(sum(tr_at >= 0) - sum(tr_at < 0)))]
tr_pl <- spplot(NY8, c("sar_trend"), at = tr_at, col = "transparent", col.regions = tr_rds, main = list(label = "Trend", cex = 0.8))
st_at <- seq(-0.16, 0.39, length.out = 21)
st_rds <- rds(sum(st_at >= 0) * 2)[ - (1:(sum(st_at >= 0) - sum(st_at < 0)))]
st_pl <- spplot(NY8, c("sar_stochastic"), at = st_at, col = "transparent", col.regions = st_rds, main = list(label = "Stochastic", cex = 0.8))
plot(tr_pl, split = c(1, 1, 2, 1), more = TRUE)
plot(st_pl, split = c(2, 1, 2, 1), more = FALSE)

#CAR MODEL
nycar <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = NY8, family = "CAR",
   listw = NYlistw)
summary(nycar)

NY8$car_trend <- nycar$fit$signal_trend
NY8$car_stochastic <- nycar$fit$signal_stochastic
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
tr_at <- seq(-1, 1.3, length.out = 21)
tr_rds <- rds(sum(tr_at >= 0) * 2)[ - (1:(sum(tr_at >= 0) - sum(tr_at < 0)))]
tr_pl <- spplot(NY8, c("car_trend"), at = tr_at, col = "transparent", col.regions = tr_rds, main = list(label = "Trend", cex = 0.8))
st_at <- seq(-0.16, 0.39, length.out = 21)
st_rds <- rds(sum(st_at >= 0) * 2)[ - (1:(sum(st_at >= 0) - sum(st_at < 0)))]
st_pl <- spplot(NY8, c("car_stochastic"), at = st_at, col = "transparent", col.regions = st_rds, main = list(label = "Stochastic", cex = 0.8))
plot(tr_pl, split = c(1, 1, 2, 1), more = TRUE)
plot(st_pl, split = c(2, 1, 2, 1), more = FALSE)


