#************#
# Distance
#***********#

#=================================#
# point coord data and data matrix
#=================================#

A <- c(40, 43)
B <- c(101, 1)
C <- c(111, 54)
D <- c(104, 65)
E <- c(60, 22)
F <- c(20, 2)

pts <- rbind(A, B, C, D, E, F)
pts

plot(pts, xlim = c(0, 120), ylim = c(0, 120), pch = 20, cex = 2,
     col = 'red', xlab = "X", ylab = "Y", las = 1)

text(pts + 5, LETTERS[1:6])


#--------------#
# dist function
#--------------#

dis <- dist(pts)
dis


#------------------------------------------#
# transform a dist matrix into a normal one
#------------------------------------------#
D <- as.matrix(dis)
round(D)


#============================#
# Distance for Lon/Lat coords
#============================#
# points were coordinates in degrees (lon/lat)

# dist() only calculates cartesian distance

# so need to use pointDistance() function

View(pointDistance)
gdis <- pointDistance(pts, lonlat = TRUE)
gdis




















