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
	library(cowplot)

### this section loads in the disparate meta-data files and concatenates them.
	### load in DrosEU data
		dat.drosEU <- read.xls("/Users/alanbergland/Documents/work/Projects/2016_DrosEU_integration/population_info/DrosEU_allYears_180607.xlsx")

		dat.drosEU.dt <- as.data.table(dat.drosEU[-1,c(2,2,5,6,7,12,13,15,16)])
		setnames(dat.drosEU.dt,
				names(dat.drosEU.dt),
				c("sampleId", "sequenceId", "country", "city", "collectionDate", "lat", "long", "season", "nFlies"))
		dat.drosEU.dt[,locality:=paste(tstrsplit(sampleId, "_")[[1]],
									tstrsplit(sampleId, "_")[[2]], sep="_")]
		dat.drosEU.dt[season=="S", season:="spring"]
		dat.drosEU.dt[season=="F", season:="fall"]
		dat.drosEU.dt[,type:="pooled"]
		dat.drosEU.dt[,collectionDate := as.POSIXct(collectionDate)]
		dat.drosEU.dt[,continent:="Europe"]
		dat.drosEU.dt[,set:="DrosEU"]


	### load in DrosRTEC data
		dat.drosRTEC <- read.xls("/Users/alanbergland/Documents/work/Projects/2016_DrosEU_integration/population_info/vcf_popinfo_Oct2018.xlsx")

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
		dat.dpgp <- read.xls("/Users/alanbergland/Documents/work/Projects/2016_DrosEU_integration/population_info/TableS2_populations.xls")

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

	### save
		save(samps, file="~/samps.Rdata")
		
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

### quality checks: are all the samples from the meta-data files in the SYNC file?
	### load header info for SYNC file
		sync.meta <- fread("/Users/alanbergland/Documents/work/Projects/2016_DrosEU_integration/population_info/DrosEU14-16-DrosRTEC-DPG.meta")
		setnames(sync.meta, names(sync.meta), c("colNum", "sequenceId", "origin"))

	### do test
		sapply(samps$sequenceId, function(x) x%in%sync.meta$sequenceId)

		samps$sequenceId[!sapply(samps$sequenceId, function(x) x%in%sync.meta$sequenceId)]
		sync.meta$sequenceId[!sapply(sync.meta$sequenceId, function(x) x%in%samps$sequenceId)]




		setkey(samps, sequenceId)
		setkey(sync.meta, sequenceId)


		samps <- merge(samps, sync.meta)



### find sites with multiple time points
	samps.ag <- samps[,list(nSamps=length(locality),
						nSpring=sum(season=="spring"),
						nFall=sum(season=="fall"),
						nTime=length(unique(collectionDate)),
						maxDelta=max(yday) - min(yday),
						lat=mean(lat),
						long=mean(long)),
				list(locality, year, continent)]

	setkey(samps.ag, locality, year)
	setkey(samps, locality, year)


### plot multi-sample populations
	multi_samp <- ggplot(data=samps[J(samps.ag[maxDelta>10])], aes(x=yday, y=lat, color=season, group=locality)) +
	geom_line() + geom_point() + facet_grid(.~year)


	single_samp <- ggplot(data=samps[J(samps.ag[maxDelta<10])], aes(x=yday, y=lat, color=season, group=locality)) +
	geom_line() + geom_point()




	ggplot(data=samps[J(samps.ag[maxDelta>20])], aes(x=yday, y=lat, color=continent, group=locality)) +
	geom_line() + geom_point() + facet_grid(.~year) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	xlab("Julian Day") + ylab("Latitude")



	samps[J(samps.ag[maxDelta>20]), popId := factor(sampleId, levels=sampleId[order(lat)])]

	ggplot(data=samps[J(samps.ag[maxDelta>20])], aes(x=yday, y=popId, color=season, group=popId)) +
	geom_line() + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	xlab("Julian Day") + ylab("Latitude")




	samps[J(samps.ag[maxDelta>20]), popId := factor(sampleId, levels=sampleId[order(lat)])]

	ggplot(data=samps, aes(x=yday, y=lat, color=season, group=sampleId)) +
	geom_line() + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	xlab("Julian Day") + ylab("Latitude")






### progression plot

	drosRTEC.plosG <- ggplot() +
			geom_line(data=samps[J(samps.ag[maxDelta>30])][set=="DrosRTEC"][year<=2011], aes(x=yday, y=lat, group=locality), color="grey") +
			geom_point(data=samps[J(samps.ag[maxDelta>30])][set=="DrosRTEC"][year<=2011], aes(x=yday, y=lat, color=season, group=locality))


	drosRTEC.plot <- ggplot() +
			geom_line(data=samps[J(samps.ag[maxDelta>30])][set=="DrosRTEC"], aes(x=yday, y=lat, group=locality), color="grey") +
			geom_point(data=samps[J(samps.ag[maxDelta>30])][set=="DrosRTEC"], aes(x=yday, y=lat, color=season, group=locality))


	dest.plot <- ggplot() +
			geom_line(data=samps[J(samps.ag[maxDelta>30])], aes(x=yday, y=lat, group=locality), color="grey") +
			geom_point(data=samps[J(samps.ag[maxDelta>30])], aes(x=yday, y=lat, color=continent, group=locality)) +
			facet_grid(year~.)

	plot_grid(drosRTEC.plosG, drosRTEC.plot, dest.plot, nrow)



### map plot

	world <- as.data.table(map_data("world"))

	samps.ag.ag <- samps.ag[,list(n=sum(nTime), lat=mean(lat), long=mean(long)), list(locality)]

	### make maps

		min.lat.eu <- 35
		max.lat.eu <- 55
		min.long.eu <- -10
		max.long.eu <- 37
		# [long>=min.long.eu & long<= max.long.eu][lat>=min.lat.eu & lat<=max.lat.eu]
		#[longitude>=min.long.eu & longitude<= max.long.eu][latitude>=min.lat.eu & latitude<=max.lat.eu]


		europe <- 	ggplot() +
					geom_polygon(data = world,
								aes(x=long, y = lat, group = group), fill="lightgrey") +
					geom_point(data = samps.ag.ag,
								aes(x=long, y=lat, size=I((n-1)/2 + 4)), alpha=.5) +
					xlab("Longitude") + ylab("Latitude") + scale_fill_manual(values="black")
		## north america
		min.lat.na <- 25
		max.lat.na <- 50
		min.long.na <- -130
		max.long.na <- -65
