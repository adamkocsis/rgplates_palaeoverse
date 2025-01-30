# Demonstration of rgplates for the Palaeoverse lecture
# 2025-01-30
# Ádám T. Kocsis

# package installation from CRAN
# install.packages("rgplates")

# suggested packages from CRAN
# install.packages(c("geojson", "httr2"))

# attachment of rgplates with sf
library(rgplates) # 0.5.0

# additional packages for the demo - installable from the CRAN!
# install.packages(("icosa", "chronosphere"))
library(icosa)
library(chronosphere)

################################################################################
# 1. Online module
########################################----------------------------------------
# A. Built-ins

# Default model: MERDITH2021 - present-day, partitioning polygons
plates0 <- reconstruct("static_polygons", age=0)
plates0

# standard sf plotting
plot(plates0$geometry, col="gray", border="black")

# actual reconstruction: Albian stage ICS 2024/12: 113.2 - 100.5 ~ 107 (106.85)
plates107 <- reconstruct("static_polygons", age=107, verbose=TRUE)
plot(plates107$geometry, col="gray", border="black")

# a different built-in feature collection
coast107 <- reconstruct("coastlines", age=107)
plot(coast107$geometry, col="gray", border=NA)

# overplot
plot(plates107$geometry, col="gray", border="white")
plot(coast107$geometry, col="black", border="white", add=TRUE)

########################################
# tested model and feature collection combinations:
data(gws)
View(gws)

# different model
coast107tc <- reconstruct("coastlines", age=107, model="TorsvikCocks2017")

# example comparison
plot(coast107$geometry, col="#AA000077", border=NA)
plot(coast107tc$geometry, col="#00AA0077", border=NA, add=TRUE)
legend("left", bty="n", legend=c("MERDITH2021", "TorsvikCocks2017"), border=NA,
	fill=c("#AA000077", "#00AA0077"))

########################################
# https://gwsdoc.gplates.org/models

# override internal checking
coast107Cao <- reconstruct("coastlines", age=107, model="CAO2024") # produces ERROR!
coast107Cao <- reconstruct("coastlines", age=107, model="CAO2024", check=FALSE)
plot(coast107Cao$geometry, col="gray90", border=NA)

########################################
# results are standard sf objects

# - Projection:
proj <- "ESRI:54009"

# plotting
coast107proj <- sf::st_transform(coast107, proj)
plot(coast107proj$geometry, col="gray90", border=NA)

# - Geometric operations
plot(coast107$geometry, col="#AA000077", border=NA)

# machine-produced data
sf_use_s2(FALSE)
coast107union <- sf::st_union(sf::st_make_valid(coast107))
plot(coast107union, col="#AA000077", border=NA)
sf_use_s2(FALSE)

########################################----------------------------------------
# B. Point reconstruction

# B1. Simple coordinates
# WGS84 longitude-latitude coordinates of the two cities
sydney<- c(151.17, -33.85)
montreal<- c(-73.61, 45.52)

# all cities in a single matrix
cities<- rbind(sydney, montreal)
cities

# Using the default MERDITH2021 model
cities107 <- reconstruct( cities, age = 107)
# returns longitude-latitude data
cities107

# using plain longlat coordinates
plot(coast107$geometry, col="gray", border=NA)
points(cities107, col="red", pch=3, cex=3, lwd=3)

########################################
# B2. More practical example: Paleobiology Database: Cretaceous Bivalves
setwd("/mnt/sky/Dropbox/Software/rgplates/doc/palaeoverse_lecture/")

# read in some example data
# https://paleobiodb.org/data1.2/occs/list.csv?datainfo&base_name=Bivalvia&interval=Cretaceous&show=coords,timecompare
pbdb <- read.csv("data/cretaceous_bivalves_2025-01-28.csv", skip=15)
str(pbdb)

# normalize collection data
# time_contain: collections datable to the precision of one GTS global stage
colls <- unique(pbdb[, c("collection_no", "lng", "lat", "time_contain")])

# stratigrahy in the time contain
unique(colls$time_contain)

# the albian subset
albian <- colls[which(colls$time_contain=="Albian"), ]
nrow(albian)

# on a present-day map
coast <- reconstruct("coastlines", age=0)

dir.create("export", showWarnings=FALSE)
png("export/albian_bivalves_present.png", width=2000, height=1000, pointsize=24)
	par(mar=c(0.1,0.1,0.1,0.1))
	plot(coast$geometry, col="gray", border=NA)
	points(albian[, c("lng", "lat")], col="red", pch=3)
dev.off()

# point reconstruction
albianRec <- reconstruct(albian[, c("lng", "lat")], age=107)
str(albianRec)

# visualize
plot(plates107$geometry, col="lightgray", border=NA)
plot(coast107$geometry, col="darkgray", border=NA, add=TRUE)
points(albianRec, col="red", pch=3)

# there are missing values here
sum(is.na(albianRec[,"paleolat"]))

# with a different model?
albianRecPM <- reconstruct(albian[, c("lng", "lat")], age=107, model="PALEOMAP")
sum(is.na(albianRecPM[,"paleolat"]))

# missing values? -> make points not inherit plate durations - DANGEROUS!
albianRec_all <- reconstruct(albian[, c("lng", "lat")], age=107, validtime=FALSE)
sum(is.na(albianRec_all[,"paleolat"]))

# add these to the map
points(albianRec_all[is.na(albianRec[,"paleolat"]), ], col="blue", pch=3, cex=2,lwd=3)


########################################
# Paleocoordinate input
plot(plates107$geometry, col="lightgray", border="white")
plot(coast107$geometry, col="gray", border="white", add=TRUE)
points(cities107, col="red", pch=3, cex=3, lwd=3)

# generate random points on the sphere
set.seed(1)
randPoints <- icosa::rpsphere(10000, output="polar")
points(randPoints, col="#00FF0033", cex=0.5, pch=16)

# spherical great circle distances from Montreal's paleocoordinates
dists <- icosa::arcdistmat(randPoints, cities107[2,, drop=FALSE])

# points closer than 4000 kms
close107 <- randPoints[dists<4000, ]
points(close107, col="#FF000088", cex=0.5, pch=16)

# where are these today?
close <- reconstruct(close107, from=107, age=0)

# replot to see the missing results
plot(plates107$geometry, col="lightgray", border="white", xlim=c(-85, 37), ylim=c(0,90))
plot(coast107$geometry, col="darkgray", border="white", add=TRUE)
points(close107, col="#FF000088", cex=0.5, pch=16)

# for which do we get missing values?
points(close107[is.na(close[,1]), ], cex=3, pch=1, col="blue")


# present-day location of the points
plot(coast$geometry, col="gray", border="white")
points(close, col="#FF000088", cex=0.5, pch=16)


########################################
# reproducibility?
reconstruct("coastlines", age=15,verbose=TRUE)

# the location of the web service
getgws()

# If GWS can run locally using docker
# setgws("http://localhost:18000/", check=FALSE)
# reconstruct("coastlines", age=15,verbose=TRUE)

################################################################################
# 2. Offline module
########################################----------------------------------------

# requires
# 1. GPlates application installation: https://www.earthbyte.org/download-gplates-2-5/
# 2. Tectonic model files

# 2A. Construct a platemodel object from files (The PaleoMAP model)
# model wrapper
pm <- platemodel(
	rotation="data/paleomap/PALEOMAP_PlateModel.rot",
	features=c(static_polygons="data/paleomap/PALEOMAP_PlatePolygons.gpml")
)

# execute reconstruction: platemodel => model
plates107PM <- reconstruct("static_polygons", age=107, model=pm)
albianRecPM <- reconstruct(albian[, c("lng", "lat")], age=107, model=pm)

# visualize
plot(plates107PM$geometry, col="gray")
points(albianRecPM, col="red", pch=3)

# 2B. Download platemodel from the Chronosphere
merdith <- chronosphere::fetch(src="GPlates", ser="MERDITH2021", ver="1.1")
merdith
# many still lack support by GPLates itself, e.g. topologies are not resolved!

# more feature collections available
merdithPoly <- reconstruct("static_polygons", age=107, model=merdith)
merdithContinents <- reconstruct("continents", age=107, model=merdith)
merdithCratons <- reconstruct("cratons", age=107, model=merdith)

# visualize
plot(merdithPoly$geometry, col="gray90", border=NA)
plot(merdithContinents$geometry, col="gray50", border=NA, add=TRUE)
plot(merdithCratons$geometry, col="black", border=NA, add=TRUE)


# Also faster!
########################################

# 3B. Multiple ages, reconstruct coordinates of Barremian, Aptian and Albian collections
# from ICS 2024/12 (Barremian bottom - Albian top)
# middle ages for three stages
bounds <- c(125.77, 121.4,  113.2, 100.5)
mids <- (bounds[2:length(bounds)-1] + bounds[2:length(bounds)]) /2
names(mids) <- c("Barremian", "Aptian", "Albian")
mids

# reconstruct the coastline feature with a model, default: enumerate=TRUE
coastMulti <- reconstruct("coastlines", age=mids, model=merdith)
str(coastMulti, 1)

# add the target age as a column
colls$map <- mids[colls$time_contain]
table(colls$map)

# get rid of data coming from not this interval
useColls <- colls[!is.na(colls$map), ]

# reconstruct - coordinate - target age pairs
usePaleo <- reconstruct(useColls[, c("lng", "lat")], age=useColls$map, model=merdith, enumerate=FALSE)
useColls <- cbind(useColls, usePaleo)
str(useColls)

pdf("export/cret_biv.pdf", width=16, height=8)
par(mar=c(0.1,0.1,2,0.1))
for(i in 1:length(mids)){
	# the current stage
	stage <- names(mids)[i]

	# make a plot
	plot(coastMulti[[i]]$geometry,
		main=paste0(stage, " ", mids[i], "Ma"),
		col="gray", border=NA)

	# visualize paleocoordinates
	points(useColls[useColls$time_contain==stage, c("paleolong", "paleolat")], pch=3, col="red")

}

dev.off()

# merge back with the rest of the occurrences
records <- merge(
	pbdb,
	useColls[, c("collection_no","paleolong", "paleolat", "map")],
	by="collection_no",
	all=TRUE)
