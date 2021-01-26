#********************#
# Spatial influence
#********************#
# an important step in spatial statistics and 
# modelling is to get a measure of spatial influence
# between geographic objs
# this can be expressed as a function of 
# adjacency or inverse distance
# and is often expressed as a spatial weights matrix


# influence is very complex and can have various ways
# to be estimated
  # e.g. the influence btw a set of polygons
  # can be expressed as having a shared border
  # or as the distance btw their centroids
  # as the lenght of a shared border, etc


#===========#
# Adjacency
#===========#
# adjacency can be thought as "touch", e.g. neighbouring countries
# can also be based on distance
# which is most common approach when analyzing point data

# define points as "adjancent" if they are within distance of 50 from each other

adj <- D < 50
adj # logical

# turn the logical adj into numeric
diag(adj) <- NA
Adj_50 <- adj * 1
Adj_50



#========================#
# Nearest neighbours
#========================#

# e.g. find the two nearest neighbours

order # arrange in an asending (default)
## 1st, for each row, re-order the cols according to the values
round(D)
cols_ord <- apply(D, 1, order)
str(cols_ord)

cols_ord_t <- t(cols_ord)

## 2nd, to get the 2 nearst neighbour
cols <- cols_ord_t[, 2:3]
cols 

# make rowcol pairs
rowcols <- cbind(rep(1:6, each = 2), as.vector(t(cols)))
head(rowcols)


## 3rd use these row-col pairs as indices to change values in matrix
#need a matrix framework with all 0 values
A_2NB <- Adj_50 * 50
A_2NB[rowcols] <- 1
A_2NB



#===============#
# weights matrix
#===============#

# adjacency express the spatial influence
# as binary value, it could also be expressed
# as a continuous value
# the simplest approach is to use inverse distance
# the further away, the lower the value

W <- 1/D
round(W, 4)

# spatial weights matrix is often row-normalized
# i.e. the sum of weights for each row in the matrix
# is the same

W[is.infinite(W)] <- NA
W

# compute row sum
r_sum <- rowSums(W, na.rm = TRUE)
r_sum

#apply(W, 1, sum, na.rm  = TRUE)

# now divid each number of W by rowsum
W <- W / r_sum

# check if each row sum up to 1
apply(W, 1, sum, na.rm = TRUE)
# yes



#==============================#
# spatial influence for polygons
#===============================#

# Above are adjacnecy for a set of points
# we look at it for polygons

library(raster)
p <- shapefile(system.file("external/lux.shp", package = "raster"))
str(p)

install.packages("sf")
install.packages("spdep")
library(sf)
library(spdep)


wr <- poly2nb(p, row.names = p$ID_2, queen = FALSE)
wr # Neighbour list
wr[1:2]

wm <- nb2mat(wr, style = "B", zero.policy = TRUE)
# the function creates an n by n weights matrix with 
# values given by the coding scheme style chosen. 
# B is the basic binary coding, 
# W is row standardised, 
# C is globally standardised, 
# while S is the variance-stabilizing coding scheme 

dim(wm)

wr[1:5]
wm[1:5, 1:11]

# compute the number of neighbours for each area / polygon
i <- rowSums(wm)
i

# in percentage of each times 
round(100 * table(i) / length(i), 1)
# 42 % of polygons have 3 neighours


# plot the links btw polygons
par(mai = c(0, 0, 0, 0))

plot(p, col = "gray", border = "blue")

xy <- coordinates(p) # retrieve
head(xy, 3) # x,y coords
dim(xy)
plot(wr, xy, add = TRUE, col = 'red', lwd = 2)


#---------------------------------------------------------#
# some alternative ways to compute the "spatial influence"
#---------------------------------------------------------#

# Distance based
wght_dist10 <- dnearneigh(xy, d1 = 0, d2 = 10)
# euclidean distance
str(wght_dist10)
wght_dist10[1:2]


wght_dist25 <- dnearneigh(xy, 0, 25, longlat = TRUE)


# Nearst neighbours
k3 <- knn2nb(knearneigh(xy, k = 3)) # retrun a neighbour list of class nb
k6 <- knn2nb(knearneigh(xy, k = 6))

k3[1:2]

# now plot all of the neighbour situations
plotit <- function(nb, lab = '') {
  plot(p, col = "gray", border = "blue")
  plot(nb, xy, col = "red", pch = 20, add = TRUE)
  text(6.3, 50.1, paste("(", lab, ")"), cex = 1.25)
}

rep(0, 4)
c(0, 0, 0, 0)
par(mfrow = c(2, 3), mai = rep(0, 4))

plotit(wr, "adjency")
plotit(wght_dist10, "weighted distance 10km")
plotit(wght_dist25, "weighted dist 25km")
plotit(k3, "3 nearest neighb")
plotit(k6, "6 nearest neighb")

