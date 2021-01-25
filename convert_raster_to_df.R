#************************************#
# Convert raster outputs to df format
#************************************#

# ref: https://gis.stackexchange.com/questions/308824/convert-raster-outputs-to-data-frame-format
# ref: https://rdrr.io/cran/raster/man/as.data.frame.html



install.packages("sp")
install.packages("raster")
library(sp)
library(raster)


r1 <- raster(matrix(1:12, 3, 4))
crs(r1) <- "+proj=longlat +datum=WGS84"

as.matrix(r1)


r2 <- r1 + 1
r3 <- r1 + 2

s <- stack(list(r1 = r1, r2 = r2, r3 = r3))
names(s)


as.data.frame(s, xy = TRUE)
as.data.frame(s, xy = TRUE, long = TRUE)
as.data.frame(s, centroids = TRUE) # no xy 

cellFromXY(s, c(0.875, 0.1666667)) # 12
cellFromXY(s, c(0.625, 0.5000000)) # 7

s[12] <- "NA"
s[12]

s[7] <- "NA"

as.data.frame(s, xy = TRUE)

as.data.frame(s, xy = TRUE, sepNA = TRUE)

