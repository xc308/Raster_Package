#**************************#
# Raster data manipulation
#**************************#
# see help("raster-package") for index of functions by topic

#=====================#
# Creating Raster* obj
#=====================#

# RasterLayer: using function raster
  # default settings: create a global raster data str
  # with a long/lat crs and 1 by 1 deg cells

  # can change these settings by providing additional arg
  # e.g. xmn, nrow, ncol, crs 
  # to change a RasterLayer to another crs, use projectRaster


# creating and changing a RasterLayer obj r 
library(raster)
x <- raster()
x

# find the resolution
res(x)

# change it 
res(x) <- 100
res(x)
x


# change ncol, which affects the resolution
ncol(x)
ncol(x) <- 18
res(x)
# [1]  22.22222 100.00000



# set the crs (define projection)
projection(x)
# [1] "+proj=longlat +datum=WGS84 +no_defs"

projection(x) <- "+proj=utm +zone=48 +datum=WGS84"
x

# above raster data only consist of geometry
# the nrows, ncols, and where the raster is located in geographic space
# no cell value associated with it

# setting and accessing values 

# start from empty raster geometry
r <- raster(ncol = 10, nrow = 10)
hasValues(r)

# assign values 
values(r) <- 1:ncell(r)
hasValues(r) # [1] TRUE


# another example
set.seed(15-01-2021)
values(r) <- runif(ncell(r))
hasValues(r) # [1] TRUE
inMemory(r)  # [1] TRUE values are in RAM 

val_r <- values(r) # a vector 1-d
values(r)[1:10]

plot(r, main = "100 cells of Raster")

# NOTE: when the ncol and/or nrow are changed
# will lose the values associated with RasterLayer if there were any
# the same applies if change the resolution directly
# but values are not lost when changing the extent (coverage)
# just like clip or enlarge the area
# so ncol or nrow not chane, but adjust the res

hasValues(r)
res(r) # [1] 36 18
dim(r) [1] 10 10  1
xmax(r) # [1] 180

# now change the maximum x coord of the extent of the RasterLayer
xmax(r) <- 0
res(r) # [1] 18 18
dim(r) # [1] 10 10  1 unchanged ncol nrow
hasValues(r) # [1] TRUE  still has values


# add the ncol, the value disappear
ncol(r) <- 5
hasValues(r) # [1] FALSE
res(r) # [1] 36 18
xmax(r) # [1] 0



# Raster package is able to deal with large file 
# as they only create a raster obj with data str
# with no values, and won't read the cell values in RAM

# if compute, data is processed in chunk

# if no filename is assigned to a function 
# or the raster output is too large, 
# the results are written to temporay file


filename <- system.file("external/test.grd", package = "raster")
r <- raster(filename)
hasValues(r) # [1] TRUE
inMemory(r) # [1] FALSE
plot(r, main = "RasterLayer from file muese")


# create multi-layer obj from RasterLayer objs
# or from files
r1 <- r2 <- r3 <- raster(nrow = 10, ncol = 10)
values(r1) <- runif(ncell(r1))
hasValues(r1)
values(r2) <- runif(ncell(r2))
hasValues(r2)
values(r3) <- runif(ncell(r3))
hasValues(r3)

# combine 3 RasterLayer obj into a RasterStack
s <- stack(r1, r2, r3)
s
nlayers(s) # [1] 3


# or combine the RasterLayer obj into a RasterBrick
b1 <- brick(r1, r2, r3)

# eqivalent to 
b2 <- brick(s)


# or create a RasterBrick from a file
filename <- system.file("external/rlogo.grd", package = "raster")
b <- brick(filename)
b

# extract a single layer from a RasterBrick or RasterStack
r_2 <- raster(b, layer = 2)


# equivalent to create from disk with band = 2
r_bd2 <- raster(filename, band = 2)
hasValues(r_bd2)

r_bd2

# =======================#
# Modifying a Raster* Obj
#========================#
# There are several functions deals with 
# modifying the spatial extent of Raster OBJ

# THE crop function lets you take a geographic subset
# of a larger raster obj

# can crop a Raster by providing an extent obj
# or another spatial obj from which an extent 
# can be extracted (objs from Raster or Sp)

# An easy wat to get an extent obj is to plot a 
# RasterLayer and then use drawExtent
# to visually determine the new exetet(bouding box)
# to provide to the crop function


# trim: crops a RasterLayer by removing the outer
# rows and cols that only contain NA values

# extend adds new rows / cols with NA values
# so could create a new RasterLayer with teh same
# extent of another, larger RasterLayer s.t.
# they can be used together in otehr functions


# merge: enable you to merge 2 or more Raster obj
# into a single new obj
# the input objs must have the same resolution and origin
# if this is not the case, need to adjust one of the raster using disaggreagte


# aggregate and disaggregate allow fo changing the resolution
# (cell size) of a Raster obj


# projectRaster(): transform values of Raster obj
# to a new obj with a different crs


r <- raster()
r[] <- 1:ncell(r)
r
r_ag <- aggregate(r, 20)
r_ag

r_disag <- disaggregate(r_ag, 20)
r_disag


r1 <- crop(r, extent(-50, 0, 0, 30))
r2 <- crop(r, extent(-10, 50, -20, 10))
m <- merge(r1, r2, filename = "test.grd", overwrite = TRUE)
plot(m)




#=========#
# Overlay
#=========#

# overlay : more efficient for large obj
  # can combine mulitple Raster obj

# mask: remove all values from one layer that are Na in another 
# cover: combines two layers by taking the values of the 1st layer
# except where these are NA


#======#
# Calc
#======#

# calc: allows you do computation for a single
# raster obj by providing a function

#==========#
# Reclasify
#==========#

# use cut or reclassify to replace ranges
# of values with single values,
# or subs to substitute single values with other values

r <- raster(ncol = 3, nrow = 2)
r[] <- 1:ncell(r)
getValues(r)

# want to set all values above 4 to NA
s <- calc(r, fun = function(x) {x[x > 4] <- NA; return(x)})
as.matrix(s) # convert a raster layer into a matrix value


# divide the 1st raster with 2 * sqrt(2nd) + 5
w <- overlay(r, s, fun = function(x, y) {x / (2 * sqrt(y)) + 5})
as.matrix(w)


# remove all values in r that are NA in w
u <- mask(r, w)
as.matrix(u)



# identify the cell values in u that are the same in s
as.matrix(s)

v <- u == s
as.matrix(v)

as.matrix(r)

# repalce NA values in w with the values in r
cvr <- cover(r, w) # combine 2 layers by taking the value in the 1st one
as.matrix(cvr) # worked!



#=========#
# Distance
#=========#

# distance: computes the shortest distance to cells that are not NA
# pointDistance: computes the shortest distance to any point in a set of points
# gridDistance: computes the distance when grid cells can be travered
# direction: computes the direction toward nearest cell that's not NA
# adjency: deteremines which cells are adjacent to other cells
# gdistance: package for advanced distance



#======================#
# saptail configuration
#======================#

# clump(): identify groups of cells that are connected
# boundaries(): identifies edges, transitions btw cell values
# area(): computes the size of each grid for unprojected rasters

r <- raster(nrow = 45, ncol = 90)
r[] <- round(runif(ncell(r)) * 3)
a <- area(r)
zonal(a, r, fun = "sum")

plot(a)



#============#
# Predictions
#============#
# two functions to make model predictions to 
# potentially very large rasters. 
  # predict: takes a multi-layer raster and a fitted model as args
    # fitted models can be various classes, e.g. lm , glm, gam, random forest
  # interpolate: for models uses coords as predictor variables, e.g. kriging, spline interpolation


#============================#
# Vector to raster conversion
#=============================#

# raster package supports point, line, polygon
# to raster conversion with rasterize function

# for vector type data (points, lines, polygons)
# objs of Spatial* classes defined in sp package are used
# but points can also be represented by a two-col matrix (x, y)

# point to raster conversion:
  # often done with the purpose to analyze the point data
  # e.g. count the number of distinct species occurs in each raster cell
  # rasterize takes a Raster obj to set the spatial extent and resolutio
  # and a function to determine how to summarize the points or attributie of each point by cell


# Polygon to raster conversion:
  # done to create a RasterLayer that can act as a mask, to set NA to a set of cells 
  # or to summarize values on a raster by zone
  # e.g. a country polygon is transferred to a raster then used to set all cells outside the country to NA
  # e.g. polygons representing administrative regions can be transferred to a raster to summarize raster values by region


# conversely, convert values of a RasterLayer to points or polygons
# using rasterToPoints and rasterToPolygons
# both functions only return values for cells that are not NA


#=======================#
# Summarizing functions
#========================#

# with the 1st arg being a Raster obj
# normal summary stats function, e.g. min, max, mean retruns a RasterLayer
# can use cellStats, if a summary of all cells of a single Raster Obj is wanted
  # freq: make a freq talbe
  # zonal: summarize a Raster obj using zones defined in a RasterLayer
  # crosstab: to cross-tabluate two RasterLayer objs

r <- raster(ncol = 36, nrow = 18)
r[] <- runif(ncell(r))
mean_all_cell <- cellStats(r, mean)


s <- r
s[] <- round(runif(ncell(s)) * 5)
s
r
zonal(r, s, 'mean') # summarize the value of r by each zone defined by s

freq(s)
#      value count
# [1,]     0    58
# [2,]     1   138
# [3,]     2   140
# [4,]     3   120
# [5,]     4   123
# [6,]     5    69

freq(s, value = 3)


ctab <- crosstab(r * 3, s)
ctab
plot(r * 3)


#=================#
# Helper functions
#=================#

# The cell number is an important concept in raster package
# Raster data can be thought of as a matrix 
# but in a RasterLayer, more commonly treated as a vector

# cells are numbered from the upper left cell
# to the upper right, then continuouing on the left side of the next row
# until the last cell at the lower-right side of the raster

# several helper functions to determine the col/ row number from a cell and vice versa
# to determine the cell number for x, y coordinates and vice versa

library(raster)

r <- raster(ncol = 36, nrow = 18)
ncol(r) # [1] 36
nrow(r) # [1] 18
ncell(r) # [1] 648


rowFromCell(r, 100) # [1] 3
# the 100th cell is on the 3rd row

colFromCell(r, 100) # [1] 28
# the 100th cell is on the 28th col

cellFromRowCol(r, 5, 5) # [1] 149
# the cell at the 5th row and 5th col is 149th cell

xyFromCell(r, 100)
#       x  y
# [1,] 95 65
# the x,y coord of the 100th cell 

cellFromXY(r, c(0, 0)) # [1] 343
# the c(0, 0) corresponds to the 343th cell

colFromX(r, 0) # [1] 19
# when x = 0 is the 19th col

rowFromY(r, 0) # [1] 10
# when y = 0 is the 10th row


#======================#
# Accessing cell values
#=======================#

# getValues(): get all values or a single row
# getValuesBlock(): read a block (rectangle) of cell values

r <- raster(system.file("external/test.grd", package = "raster"))
r_50 <- getValues(r, 50) # get the 50th row value

ncol(r) # [1] 80
nrow(r) # [1] 115

r_50 # a vector
r_50[35: 39]
# [1] 743.8288 706.2302 646.0078 686.7291 758.0649

r
getValuesBlock(r, 50, 1, 35, 5)
# getValuesBlock(x, row=1, nrows=1, col=1, ncols=(ncol(x)-col+1), format='', ...)
#[1] 743.8288 706.2302 646.0078 686.7291 758.0649
# the block starts from the 50th row, lasts for 1 row, 
  # and starts from the 35th col, lasts for 5 columns
# which is the same as :
  # r_50 <- getValues(r, 50)
  # r_50[35: 39]



# also read values using cell numbers or coords
# using extract method
# but 1st need to get the cell numbers or coords
# then extract the corresponding values


cell_nums <- cellFromRowCol(r, 50, 35:39)
extract(r, cell_nums)
# [1] 743.8288 706.2302 646.0078 686.7291 758.0649

xy_coords <- xyFromCell(r, cell = cell_nums)
extract(r, xy_coords)
# [1] 743.8288 706.2302 646.0078 686.7291 758.0649


# can also extract using SpatialPolygons* or SpatialLines*
# Default approach for extracting raster values with polygons 
# is that a polygon has to cover the center of a cell, for that cell to be included
# but can use argument "weights=TRUE" 
# get both the cell values,and the % of each cell that's covered by the polygon
# so can apply e.g. a 50% area covered threshold 
# or compute an area-weighted average

# for lines, any cell that's crossed by a line is inclueded
# for lines and points, a cell is included when it is below or to the right of the line/point

# can also use standard R indexing to access values
# to replace values, assign new values in a raster obj

r[cell_nums]
# [1] 743.8288 706.2302 646.0078 686.7291 758.0649

r[1:4] # [1] NA NA NA NA


# setting new raster values of a file, 
# then the connection to that file will be lsot
# as it is now different from the original file


filename(r)
# [1] "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/raster/external/test.grd"

r[2:3] <- 10
r[1:4] # [1] NA 10 10 NA

filename(r) # [1] ""


# above, values are retrieved using cell numbers
# that's a raster is represented using 1-d vector

# values can also be inspected using 2-d matrix
# 

r[1, ] # the 1st row of raster r
r[50, ]
r[, 2]

r[1:3, 1:3]

r[1:3, 1:3, drop = FALSE] # retain matrix str

# but accessing values through this type of indexing shoule
# be avoided inside functions as less efficient than 
# via getValues()



























































































































######################################
assign_val <- function(r = r) {
  values(r) <- runif(ncell(r))
}

lapply(list(r1, r2, r3), assign_val)
######################################











































