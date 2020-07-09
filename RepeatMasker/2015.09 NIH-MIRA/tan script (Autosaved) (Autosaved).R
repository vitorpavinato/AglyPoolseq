### libraries
	library(data.table)
	library(foreach)
	library(lme4)
	
	
### data
	load("/mnt/nescent/2014-mel-seasonality/analysis/analyses_data/00d-sfw_data_all.RData")
	setkey(sfw_data_all, chrom, pos)
	
	
	plotFun <- function(r) {
		#tan <- sfw_data_all[J(pop.pig[rank==r]$variant_id)]
		tan <- sfw_data_all[J(r)]
		
		tan$pop_factor <- factor(tan$pop_name, levels=c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12", "BA_12", "VI_12"))
		tan$season_factor <- factor(tan$season, levels=c("spring", "fall"))
		
		tan$x <- as.numeric(tan$pop_factor) + as.numeric(tan$season_factor)/3	
		tan$se <- sqrt(tan$af*(1-tan$af)/tan$dp)
		
		tan <- na.omit(tan)
		
		tan.ag <- tan[,list(delta=abs(af[season_factor=="spring"] - af[season_factor=="fall"])), list(pop_factor=pop_factor)]
		
		
		tan.glm <- summary(glmer(af~season_factor + (1|pop_factor), tan, weights=tan$dp, family=binomial()))
		
		setkey(tan, pop_factor, season_factor)
	
		
		tan$exp_af <- tan$af
		tan[season_factor=="fall"]$exp_af <- plogis(qlogis(tan[season_factor=="spring"]$af) + coef(tan.glm)[2,1])
		
		
	### plot
		plot(af~x, tan, type="n", ylim=c(min(tan$af-2*tan$se,1), max(tan$af-2*tan$se,1)), main=paste(tan$chrom[1], tan$pos[1], sep=" "))
		
		for(i in c("VA_12", "PA_09", "PA_10", "PA_11", "PA_12", "MA_12", "NY_12", "WI_12", "BA_12", "VI_12")) {
			
			arrows(x0=na.omit(tan[pop_factor==i]$x), x1=na.omit(tan[pop_factor==i]$x), y0=na.omit(tan[pop_factor==i]$af+1.96*tan[pop_factor==i]$se),
					na.omit(tan[pop_factor==i]$af-1.96*tan[pop_factor==i]$se), code=3, length=.01, angle=90)
	
			
			points(af~x, na.omit(tan[pop_factor==i]), type="b", pch=19, col=c("blue", "red"), cex=1.5)
			points(exp_af~x, na.omit(tan[pop_factor==i]), type="b", pch=21, col=c("blue", "red"), cex=1.5, lty="dashed")
			
		}
	}
	
	pdf("~/out.pdf", h=5, w=5)
	for(i in 1:500) try(plotFun(i), silent=T)
	dev.off()	
		
		
		