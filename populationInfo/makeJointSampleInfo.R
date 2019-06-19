### Make joint popInfo file for DEST
### Final data has columsn: sampleId, country, city, collectionDate, lat, long, season, nFlies, locality, type (inbred/pooled), continent
### Alan Bergland, Oct 3, 2018

### libraries
	library(data.table)
	library(gdata)
	library(cowplot)
	library(data.table)
	library(foreach)
	library(ggplot2)
	library(ggmap)
	library(maps)
	library(mapdata)

### this section loads in the disparate meta-data files and concatenates them.
	### load in DrosEU data
		dat.drosEU <- read.xls("./populationInfo/DrosEU_allYears_180607.xlsx")

		dat.drosEU.dt <- as.data.table(dat.drosEU[-1,c(2,2,5,6,7,12,13,15,16)])
		setnames(dat.drosEU.dt,
				names(dat.drosEU.dt),
				c("sampleId", "sequenceId", "country", "city", "collectionDate", "lat", "long", "season", "nAutosomes"))
		dat.drosEU.dt[,locality:=paste(tstrsplit(sampleId, "_")[[1]],
									tstrsplit(sampleId, "_")[[2]], sep="_")]
		dat.drosEU.dt[season=="S", season:="spring"]
		dat.drosEU.dt[season=="F", season:="fall"]
		dat.drosEU.dt[,type:="pooled"]
		dat.drosEU.dt[,collectionDate := as.POSIXct(collectionDate)]
		dat.drosEU.dt[,continent:="Europe"]
		dat.drosEU.dt[,set:="DrosEU"]
		dat.drosEU.dt[,nFlies:=nAutosomes/2]

	### load in DrosRTEC data
		dat.drosRTEC <- read.xls("./populationInfo/vcf_popinfo_Oct2018.xlsx")

		dat.drosRTEC.dt <- as.data.table(dat.drosRTEC[,c(1, 4, 9, 7, 12, 10, 11, 6, 16, 3)])
		setnames(dat.drosRTEC.dt,
				names(dat.drosRTEC.dt),
				c("sampleId", "sequenceId", "country", "city", "collectionDate", "lat", "long", "season", "nFlies", "locality"))
		dat.drosRTEC.dt[,type:="pooled"]
		dat.drosRTEC.dt[,collectionDate := as.POSIXct(collectionDate)]
		dat.drosRTEC.dt[long>0,continent:="Europe"]
		dat.drosRTEC.dt[long<0,continent:="NorthAmerica"]
		dat.drosRTEC.dt[,set:="DrosRTEC"]

	### load in DPGP data
		###http://johnpool.net/TableS2_populations.xls
		dat.dpgp <- read.xls("./populationInfo/TableS2_populations.xls")

		dat.dpgp.dt <- as.data.table(dat.dpgp[-c(1:4),c(1,1, 2,3,4,6,7,9)])
		setnames(dat.dpgp.dt,
				names(dat.dpgp.dt),
				c("sampleId", "sequenceId", "country", "city", "collectionDate", "lat", "long", "nFlies"))

		dat.dpgp.dt <- dat.dpgp.dt[sampleId%in%c("CO", "GA", "GU", "NG", "ZI")]
		dat.dpgp.dt[,collectionDate := paste(tstrsplit(collectionDate, "/")[[2]], tstrsplit(collectionDate, "/")[[1]], 15, sep="/")]
		dat.dpgp.dt[,collectionDate := as.POSIXct(collectionDate)]

		dat.dpgp.dt[sequenceId=="ZI", sequenceId:="dpgp3"]
		dat.dpgp.dt[,season:=NA]
		dat.dpgp.dt[,locality:=sampleId]
		dat.dpgp.dt[,type:="inbred"]
		dat.dpgp.dt[,continent:="Africa"]
		dat.dpgp.dt[,set:="dpgp"]

	### combine
		samps <- rbind(rbind(dat.drosEU.dt, dat.drosRTEC.dt), dat.dpgp.dt)
		samps[,lat := as.numeric(as.character(lat))]
		samps[,long := as.numeric(as.character(long))]
		samps[,year := year(collectionDate)]
		samps[,yday := yday(collectionDate)]


	### get GHCND site
		gh <- fread("~/ghcnd-stations.csv", header=F, fill=T)
		setnames(gh, names(gh), c("stationName", "lat", "long"))

		o <- foreach(i=1:dim(samps)[1], .combine="rbind")%do%{
			d <- sqrt((samps[i]$lat - gh$lat)^2 + (samps[i]$long - gh$long)^2)
			o <- gh[which.min(d)]
			o[,dist:=min(d)]
			o[,sampleId:=samps[i]$sampleId]
		}
		setkey(o, sampleId)
		setkey(samps, sampleId)

		samps <- merge(samps, o[,-c("lat", "long"), with=F])

	### save
		write.csv(samps, "./populationInfo/samps.csv", quote=F, row.names=F)
