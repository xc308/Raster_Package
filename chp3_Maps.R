#******#
# Maps
#******#

library(raster)
b <- brick(system.file("external/rlogo.grd", package = "raster"))
p <- shapefile(system.file("external/lux.shp", package = "raster"))
p
# class       : SpatialPolygonsDataFrame 
# features    : 12 
# extent      : 5.74414, 6.528252, 49.44781, 50.18162  (xmin, xmax, ymin, ymax)
# crs         : +proj=longlat +datum=WGS84 +no_defs 
# variables   : 5
# names       : ID_1,     NAME_1, ID_2,   NAME_2, AREA 
# min values  :    1,   Diekirch,    1, Capellen,   76 
# max values  :    3, Luxembourg,   12,    Wiltz,  312 


r <- raster(p, res = 0.01)
r
values(r) <- 1:ncell(r)
hasValues(r)

# set the those are NA in Polygons to NA in r
r <- mask(r, p)
# now raster layer r has geo shape

# if use plot to plot rasterlayer, 
# it calls the function 'rasterImage' (add legend by default)
# using code from 
View(fields::image.plot)


#---------------------------------------------------------#
install.packages("fields")
library(fields)
image.plot()

# horizontal	
# If false (default) legend will be a vertical strip on the right side. 
# If true the legend strip will be along the bottom.

# Note that you can use image to draw a bunch of images 
# and then follow with image.plot and legend.only=TRUE to add a common legend. 
# (See examples below.)
#--------------------------------------------------------------------------------#


# After plotting a RasterLayer, can add vector spatial data
# (points, lines, polygons)
# can use basic functions points, lines, polygons
# or plot(obj, add = TRUE) if using Spatial* obj defined in sp

# when plot multi-layer Raster* obj
# all layers are plotted up to 16

plot(r)
plot(p, add = TRUE)
plot(b)



levels()

quantile()

ppoint
order

pretty









