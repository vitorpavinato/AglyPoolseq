setkey(pOut, scaffold, pos)
setkey(fq, scaffold, pos)


fsq <- fq[,names(x)[grepl(".alt", names(x))], with=F] / fq[,names(x)[grepl(".rd", names(x))], with=F]
pOut$lowFq <- apply(fsq, 1, FUN=function(x) any(x<=0.01 | x>=.99))

pOut$rank[pOut$use & !pOut$lowFq] <- rank(pOut$p.glm.new[pOut$use & !pOut$lowFq])

x <- fq[pOut$rank<=1]
x.fq <- x[,names(x)[grepl(".alt", names(x))], with=F]/ x[,names(x)[grepl(".rd", names(x))], with=F]
x.fq <- cbind(x.fq[,names(x.fq)%in%paste(locales[locales$trt_noSoF=="fish"]$pop, ".alt", sep=""),with=F],	
	x.fq[,names(x.fq)%in%paste(locales[locales$trt_noSoF=="midge"]$pop, ".alt", sep=""),with=F])




plot(1~1, type="n", ylim=c(0,1), xlim=c(1,8))
for(i in 1:100) {
	temp <- as.matrix(x.fq[i])
	
	if(mean(temp[1:4], na.rm=T)<mean(temp[5:8], na.rm=T)) {
		points(as.numeric(as.matrix(temp))~c(1:8), type="l", lwd=.75, col="lightgrey")
	} else {
		points(I(1-as.numeric(as.matrix(temp)))~c(1:8), type="l", lwd=.75, col="lightgrey")
	}
}
	