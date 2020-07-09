#### libraries
	library(data.table)
	library(foreach)
	
#### data
	load("/mnt/nescent/2014-mel-seasonality/analysis/analyses_data/00d-sfw_data_all.RData")
	setkey(sfw_data_all, variant_id)
	
	load("/mnt/nescent/2014-mel-seasonality/analysis/analyses_data/00d-variant_ids.RData")
	setkey(variant_ids, chrom, pos)
	
	sfw.ag <- sfw_data_all[,list(fhat=mean(af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12")]),
								fhat_pa=mean(af[pop_name%in%c("PA_09", "PA_10", "PA_11")]),
								fvar=var(af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12")])),
							list(variant_id)]
	
	
	
	## chill coma
		chill <- read.delim("/mnt/Alan/new_alignments/vcf2/flyland/chillcoma.allele.cov", header=F, skip=1, sep=" ", as.is=T)
		names(chill) <- c("chrom", "pos", "ref", "major", "minor", "slow_major", "slow_minor", "fast_major", "fast_minor")
		chill$slow_ref_fq <- chill$slow_major/(chill$slow_major + chill$slow_minor)
		chill$slow_ref_fq[chill$ref==chill$minor] <- 1 - chill$slow_ref_fq[chill$ref==chill$minor]
		chill$slow_rd <- (chill$slow_major + chill$slow_minor)
	
		chill$fast_ref_fq <- chill$fast_major/(chill$fast_major + chill$fast_minor)
		chill$fast_ref_fq[chill$ref==chill$minor] <- 1 - chill$fast_ref_fq[chill$ref==chill$minor]
		chill$fast_rd <- (chill$fast_major + chill$fast_minor)
	
		chill$z <- (chill$fast_ref_fq - chill$slow_ref_fq) / sqrt((chill$fast_ref_fq + chill$slow_ref_fq)/2 * 
					(1-(chill$fast_ref_fq + chill$slow_ref_fq)/2) * (1/600 + 1/chill$slow_rd + 1/chill$fast_rd))
		chill$p <- 2*(pnorm(abs(chill$z), 0, 1, lower.tail=F))
	
		chill <- as.data.table(chill)
		setkey(chill, chrom, pos)
	
		chill.6d <- chill[variant_ids]
		chill.6d$rank <- rank(chill.6d$p)
	
		
	## starv resistance
		starv <- read.delim("/mnt/Alan/new_alignments/vcf2/flyland/starvation.allele.cov", header=F, skip=1, sep=" ", as.is=T)
		names(starv) <- c("chrom", "pos", "ref", "major", "minor", "control_major", "control_minor", "starved_major", "starved_minor")
		starv$control_ref_fq <- starv$control_major/(starv$control_major + starv$control_minor)
		starv$control_ref_fq[starv$ref==starv$minor] <- 1 - starv$control_ref_fq[starv$ref==starv$minor]
		starv$control_rd <- (starv$control_major + starv$control_minor)
	
		starv$starved_ref_fq <- starv$starved_major/(starv$starved_major + starv$starved_minor)
		starv$starved_ref_fq[starv$ref==starv$minor] <- 1 - starv$starved_ref_fq[starv$ref==starv$minor]
		starv$starved_rd <- (starv$starved_major + starv$starved_minor)
	
		starv$z <- (starv$starved_ref_fq - starv$control_ref_fq) / sqrt((starv$starved_ref_fq + starv$control_ref_fq)/2 * 
					(1-(starv$starved_ref_fq + starv$control_ref_fq)/2) * (1/600 + 1/starv$control_rd + 1/starv$starved_rd))
		starv$p <- 2*(pnorm(abs(starv$z), 0, 1, lower.tail=F))
	
		starv <- as.data.table(starv)
		setkey(starv, chrom, pos)
	
		starv.6d <- starv[variant_ids]

		starv.6d$rank <- rank(starv.6d$p)

#### functions

	Fst <- function(spring, fall) {
		fhat <- spring/2 + fall/2
		Htot <- 2*fhat*(1-fhat)
		
		Hwith <- spring*(1-spring) + fall*(1-fall)
		
		max((Htot-Hwith)/Htot)
	}

	matchControl <- function(targetId, param, paramId=NULL, exclude=NULL, nBoot=5, method=2) {
			if(method==1) {
				len <- length(targetId)
				if(nBoot>0) {
					out <- cbind(targetId, do.call("rbind", foreach(i=1:len, .export=c("param", "targetId", "nBoot", "len", "exclude"))%dopar% {
						print(paste(i, " / ", len, sep=""))
				#			potentialIds <- param[J(param[param$id==targetId[i], -dim(param)[2], with=FALSE])]$id  ## try to make this line vectorized
						potentialIds <- param[J(param[id==targetId[i], -dim(param)[2], with=FALSE]), nomatch=0]$id  ## try to make this line vectorized

						potentialIds <- potentialIds[!potentialIds%in%targetId]
	
						if(!is.null(exclude)) {
							potentialIds <- potentialIds[!exclude[potentialIds]]
						}
	
						ret <- NULL
						if(length(potentialIds)>0) {
							ret <- as.numeric(sample(as.character(potentialIds), nBoot, replace=TRUE))
						} else {
							ret <- rep(NA, nBoot)
						}
						ret
					}))
				} else {
					out <- matrix(targetId, ncol=1)
				}
			} else if(method==2) {
				origParamKey <- key(param)
				setkey(param, id)
				targetDt <- param[J(targetId)]
		
				setkeyv(param, origParamKey)
				setkeyv(targetDt, origParamKey)
		
				setnames(param, "id", "id.control")
				setnames(targetDt, "id", "id.target")
			
				ret <- param[targetDt, nomatch=NA, allow.cartesian=T]
			
			
				ret <- ret[!ret$id.control%in%targetId]
				ret$id[is.na(ret$id.control)]<-0
		
				sampRet <- ret[,as.numeric(sample(as.character(id.control), nBoot, replace=T)),id.target]
			
				ret <- matrix(sampRet$V1, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)
		
		
				out <- cbind(matrix(sampRet$id.target, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)[,1], ret)			
				ret[ret==0] <- NA
		
			} else if(method==3) {
				targetDt <- paramId[targetId]
		
				setkeyv(targetDt, key(param))
		
				ret <- param[J(targetDt), nomatch=NA, allow.cartesian=TRUE]
				#ret <- param[J(targetDt), nomatch=NA]
		
		
				ret <- ret[!ret$id%in%targetId]
		
				if(!is.null(exclude)) {
					ret <- ret[!ret$id%in%exclude]		
				}
		
		
				if(dim(ret)[1]>0) {
					ret$id[is.na(ret$id)]<-0
					sampRet <- ret[,as.numeric(sample(id,nBoot, replace=T)),id.1]
		
					ret <- matrix(sampRet$V1, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)
		
		
					out <- cbind(matrix(sampRet$id.1, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)[,1], ret)
					out[out==0] <- NA
				} else {
					out <- matrix(c(targetId, rep(NA, nBoot)), nrow=1)
				}
			}
			
			out
		}

	st <- function(spring, fall) {
		all(sum(sign(spring-fall))==length(spring) | sum(sign(spring-fall))==-length(spring))
	}

### analysis
	
	fst.control.mat <- matchControl(targetId=chill.6d[rank<50]$variant_id,
									param=data.table(fhat=round(sfw.ag$fhat_pa, 3),
													fvar=round(sfw.ag$fvar, 1),
													id=sfw.ag$variant_id,
													key="fhat"),
									nBoot=100)
	
	o <- foreach(i=1:dim(fst.control.mat)[2])%dopar%{
	
		sfw_data_all[J(fst.control.mat[,i]),
						list(Fst=Fst(spring=af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12") & season=="spring"],
									fall=af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12") & season=="fall"]),
							fhat=mean(af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12")]),
							st = st(spring=af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12") & season=="spring"],
									fall=af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12") & season=="fall"]),
							boot=i),
						list(variant_id)]
	}

	o <- rbindlist(o)
	plotmeans(st~boot, o)
	
	
	
	
	o.ag <- do.call("rbind", foreach(i=2:dim(fst.control.mat)[2], .errorhandling="remove")%do%{
		data.table(fhat=mean(o[boot==1]$fhat / o[boot==i]$fhat), fst=mean(o[boot==1]$Fst - o[boot==i]$Fst))
	})
	
	hist(o.ag$fst)
	hist(o.ag$fhat)


	
	control <- sfw_data_all[J(starv.6d[sample(500000, dim(target)[1])]$variant_id),
					list(Fst=Fst(spring=af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12") & season=="spring"],
								fall=af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12") & season=="fall"]),
								fhat=mean(af[pop_name%in%c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12")]),
								class="control"),
					list(variant_id)]
	control <- control[order(Fst)]
	
	plot(target$Fst~control$Fst)
	abline(0,1, col="red")
	
	
	t.test((target$Fst), (control$Fst))
	t.test(target$fhat, control$fhat)
























