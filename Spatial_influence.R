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
































