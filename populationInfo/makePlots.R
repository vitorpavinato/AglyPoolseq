### Make basic plots for DEST samples
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


### load data
  samps <- fread("./populationInfo/samps.csv")

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

	multi_sample <- ggplot() +
	geom_line(data= samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, linetype=continent)) +
	geom_point(data=samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, color=season)) +
	facet_grid(.~year) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="vertical") +
	scale_x_date(date_labels = "%b", limits = as.Date(c(110,355), origin = as.Date("2018-01-01"))) +
	xlab("Collection Date") + ylab("Latitude")

	ggsave(multi_sample, file="~/multiSample.pdf")


### num. flies plot

	samps[,nFlies:=as.numeric(as.character(nFlies))]
	samps[,sampleId:=as.character(sampleId)]
	samps[,sampleId:=factor(sampleId,
							levels=c(samps[set=="DrosEU"]$sampleId,
						  			 samps[set=="DrosRTEC"][order(nFlies)]$sampleId,
								 	 samps[set=="dpgp"]$sampleId))]


	nFlies.plot <- ggplot(data=samps, aes(x=sampleId, y=nFlies, color=set)) + geom_point() +
					theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +
					ylab("Num. flies sampled")

	ggsave(nFlies.plot, file="~/numFlies.pdf")


