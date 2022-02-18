library(maps)
library(mapdata)
library(maptools)
library(plotrix)

## Load the geo locations of each site:
gps <- read.csv(file = "filedPools_aggregated_geocoord_biotypes.csv", header = F, sep = ",")
colnames(gps) <- c("loc", "lat", "lon")
gps <- gps[-c(2,4,6,8,10,12),]

# MAP
pdf("map.pdf",         # File name
    width = 11, height = 8.50, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
)
map("state", interior = FALSE, xlim=c(-105,-70), ylim=c(38,52))
map("state", boundary = FALSE, lty = 1, add = TRUE)
map.axes(cex.axis=0.9)
map.scale(x=-104, y=51, relwidth = 0.15, metric = TRUE, ratio = TRUE, cex=1)
points(gps$lon, gps$lat, col="black", pch=19)
dev.off()





