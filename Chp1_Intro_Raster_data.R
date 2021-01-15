# ************#
# Raster Data 
#*************#

#=============#
# RasterLayer
#=============#

# create a RasterLayer

# raster require sp
install.packages("sp")
install.packages("raster")
library(raster)
library(sp)

r <- raster(ncol = 10, nrow = 10, xmn = -150, xmx = -80,
       ymn = 20, ymx = 60)
r
# obj r has a skeleton of a raster data set
# only knows its location, resolution, 
# but has no values associated with it


# assign some values
# a vector with length equal to the number of cells of RasterLayer
values(r) <- runif(ncell(r))
r

# can plot this RasterLayer obj
par(mai = c(0.3, 0.3, 0.3, 0.3))
plot(r)

# add polygon and points
lon <- c(-116.8, -114.2, -112.9, -111.9, -114.2, -115.4, -117.7)
lat <- c(41.3, 42.9, 42.4, 39.8, 37.6, 38.3, 37.6)
lonlat <- cbind(lon, lat)

polys <- spPolygons(lonlat, crs = '+=proj=longlat +datum=WGS84')
# Create spatialPolygons* obj

points(lonlat, col = 'red', cex = 3, pch = 20)
plot(polys, border = 'blue', lwd = 2, add = TRUE)


#=============================#
# RasterStack and RasterBrick
#==============================#
# make a RasterStack from multiple layers

r2 <- r * r
r3 <- sqrt(r)

s <- stack(r, r2, r3)
s

plot(s)

# make a RasterBrick from a RasterStack
b <- brick(s)
b


#=====================#
# Reading Raster files
#=====================#

f <- system.file("external/rlogo.grd", package = "raster")
f

# read
r1 <- raster(f)
r1 # the 1 out of 3 bands

# request other bands
r2 <- raster(f, band = 2)
r2

# want all layers in a single obj
# can use brick
b <- brick(f)
b
# class      : RasterBrick 
# dimensions : 77, 101, 7777, 3  (nrow, ncol, ncell, nlayers)

# these approach holds for other raster file formats
# e.g. GeoTiff, NetCDF, Imagine, ESRI Grid



#=====================#
# Writing Raster file
#=====================#

# use writeRaster to write raster data
# must provide a Raster* obj and a file name
# the file format will be guessed from the filename extension
# if it doesn't work, cna provide arg format = GTIFF



# require Rgdal

install.packages("rgdal")
library(rgdal)

# origianlly, s is a stack of r, r1, r2
x <- writeRaster(s, "output.tif", overwrite = TRUE)
# write the Raster* obj to a file
# x now is a RasterBrick
x
























































