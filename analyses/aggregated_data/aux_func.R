#####
### AUXILIARY FUNCTIONS
####

###
###
### --- ALLELE FREQUENCY AND DIVERSITY STATISTISC ---
###
###

#' COMPUTE MAXIMUM LIKELIHOOD IMPUTED SAMPLE ALLELE COUNT
#'
#' This function outputs the reference allele count with
#' a maximum likelihood imputation approach.
#' @param x The the poolseq data from poolfstat (class pooldata)
#' @return a vector of reference allele counts for each pool in x
#' @export

imputedRefMLCount <- function(x)
{
  coverage = x@readcoverage
  counts = x@refallele.readcount
  
  nbr.gene.copies = x@poolsizes[1]
  nbr.loci = x@nsnp
  nbr.pops = x@npools
  
  # Receive the outputs
  refcount.matrix = matrix(nrow=nbr.loci, ncol = nbr.pops)
  altcount.matrix = matrix(nrow=nbr.loci, ncol = nbr.pops)
  haplocount.matrix = matrix(nrow=nbr.loci, ncol = nbr.pops)
  
  for (k in 1:nbr.pops)
  {
    reads.allele.1 <- counts[,k]
    reads.allele.2 <- coverage[,k] - counts[,k]
    total.counts <- coverage[,k]
    
    allele.counts.1 <- vector(length = nbr.loci, mode = "integer")
    allele.counts.2 <- vector(length = nbr.loci, mode = "integer")
    haploid.counts  <- vector(length = nbr.loci, mode = "integer")
    
    likelihood <- matrix(nrow = nbr.loci, ncol = (nbr.gene.copies + 1))
    
    for (i in 1:nbr.loci) {
      if (total.counts[i] == 0){
        allele.counts.1[i] <- NA
        allele.counts.2[i] <- NA
        haploid.counts[i]  <- NA
        
      } else if ( total.counts[i] <= nbr.gene.copies) {
        allele.counts.1[i] <- reads.allele.1[i]
        allele.counts.2[i] <- reads.allele.2[i]
        haploid.counts[i] <- total.counts[i]  
        
      } else {
        if (reads.allele.1[i] == 0 | reads.allele.1[i] == total.counts[i]){
          allele.counts.1[i] <- reads.allele.1[i]
          allele.counts.2[i] <- reads.allele.2[i]
          haploid.counts[i] <- total.counts[i]  
          
        } else {
          for (j in 1:(nbr.gene.copies + 1)) 
          {
            likelihood[i,j] <- dbinom(x = reads.allele.1[i], size = total.counts[i],
                                      prob = (j - 1) / nbr.gene.copies, log = FALSE)
          } # end of for i,j
          
          allele.counts.1[i] <- (which(likelihood[i,] == max(likelihood[i,])) - 1)
          allele.counts.2[i] <- (nbr.gene.copies - allele.counts.1[i])
          haploid.counts[i] <- nbr.gene.copies
        } # end of inside if
          
      } # end of if
      
    } # end of for i
    
    refcount.matrix[,k] <- allele.counts.1
    altcount.matrix[,k] <- allele.counts.2
    haplocount.matrix[,k] <- haploid.counts
    
  } # end of FOR k
  
  out_list <- list(ref_count = refcount.matrix, 
                   alt_count = altcount.matrix,
                   hap_count = haplocount.matrix)
  
  return(out_list)
  
}# end of function definition

## COMPUTE GENOME-WIDE MEAN THETA FOR EACH SAMPLE
thetaHE <- function(x)
{
  l <- sum(!is.na(x))
  return((x/(1-x))/l)
}

## COMPUTE OVERALL HETEROZYGOSITY ACROSS SAMPLES - TOTAL HETEROZYGOSITY H_T
#'
#' For each SNP, it computes the heterozygosity across samples (vector)
#' @param x A vector with reference allele frequency of each sample
#' @return A scalar 
#' @export

totalHE <- function(x)
{
  p <- mean(as.numeric(x), na.rm = T)
  return(2*p*(1 - p))
}

## COMPUTE AVERAGE WITHIN SAMPLE LOCUS GENETIC DIVERSITY (h_0)
#'
#' For each SNP, It calculates the within locus pi
#' @param x A vector with reference allele frequency of each sample
#' @param Nindiv The the number of individuals in the pools
#' @return DataFrame with 1 - F_0 (within locus Pi) and F_0
#' @export

avgWithinLocusPi <- function(x, Nindv = 5)
{
  M = Nindv
  r = length(x[!is.na(x)])
  X = sum(x^2, na.rm = T) + sum((1-x)^2, na.rm = T)
  F_0 = ((2 * M * X) - r)/(((2*M)-1) * r)
  
  res = data.frame(pi_within = (1 - F_0), F_0 = F_0)
  return(res)  
    
}

## COMPUTE AVERAGE BETWEEN SAMPLE LOCUS GENETIC DIVERSITY (H_B OR H_1)
#'
#' For each SNP, It calculates the between locus pi
#' @param x A vector with reference allele frequency of each sample
#' @param Nindiv The the number of individuals in the pool
#' @return DataFrame with 1 - F_1 (between locus Pi) and F_1
#' @export

avgBetweenLocusPi <- function(x, Nindv = 5)
{
  M = Nindv
  r = length(x[!is.na(x)])
  X = sum(x^2, na.rm = T) + sum((1-x)^2, na.rm = T)
  Y = sum(x, na.rm = T)^2 + sum(1-x, na.rm = T)^2
  F_1 = (Y - X)/(r * (r - 1))
  
  res = data.frame(pi_between = (1 - F_1), F_1 = F_1)
  return(res)  
  
}

## COMPUTE GENETIC DIVERSITY 
#' from Smadja et al 2012
#' Compute for all SNPs in the SNP matrix: 
#' within-sample heterozygosity (h_0),
#' between-sample heterozygosity (H_B or H_1),
#' HE across samples (H_T), 
#' the FST
#' @param x a matrix with all ML reference allele frequencies for all samples/loci
#' @return a list
#' @export

geneticDiversity <- function(x, Nindv=5)
{
  snpInfo = x[, c(1:4)]
  x = x[, -c(1:4)]
  
  t_0 <- do.call(rbind, apply(x, 1, function(x) avgWithinLocusPi(x, Nindv = Nindv)))
  t_1 <- do.call(rbind, apply(x, 1, function(x) avgBetweenLocusPi(x, Nindv = Nindv)))
  h_t <- apply(x, 1, function(x) totalHE(x))
  
  beta = (t_0$F_0 - t_1$F_1)/(1 - t_1$F_1)
  
  return(list(snpInfo = snpInfo,
              F_0 = t_0$F_0,
              F_1 = t_1$F_1,
              pi_within = t_0$pi_within,
              pi_between = t_1$pi_between,
              H_T = h_t,
              fst = beta
              ))
  
}

## COMPUTE WITHIN-SAMPLE GENETIC DIVERSITY (pi_within)
#'
#' It calculates the within population pi for each locus 
#' @param x The scalar or vector containing the ML imputed reference allele frequencie
#' @param pool_size The 2x the number of individuals in each pool
#' @return The scalar or vector of pi_within
#' @export

withinLocusPi <- function(x, pool_size = 10)
{
  r = x
  a = 1-x
  pi <- (pool_size/(pool_size - 1)) * (r * a * 2)
  return(pi)
}

## COMPUTE F_ST
#' from Nouhaud et al 2018 Identifying genomic hotspots of differentiation and candidate
#' genes involved in the adaptive divergence of pea aphid host races. Molecular Ecology.
#' @param DATA_SIZE a matrix with all pool sizes for all samples/loci
#' @param DATA_FREQ a matrix with all ML reference allele frequencies for all samples/loci
#' @return a list
#' @export

FST_WC <- function(DATA_SIZE,DATA_FREQ){
  Nrace=dim(DATA_SIZE)[2] ; SumNi=rowSums(DATA_SIZE)
  Nic=DATA_SIZE-(DATA_SIZE**2)/SumNi ; Nc=rowSums(Nic)/(Nrace-1)
  
  MSG=(rowSums(DATA_FREQ*(1-DATA_FREQ)*DATA_SIZE)) /(SumNi-1)
  PA=rowSums(DATA_FREQ*DATA_SIZE)/SumNi
  MSP=(rowSums(DATA_SIZE*((DATA_FREQ-PA)**2)))/(Nrace-1)
  NUMERATOR=(MSP-MSG) ; DENOMINATOR=(MSP+(Nc-1)*MSG)
  
  FST=NUMERATOR/DENOMINATOR ; #FST[FST<0]=0
  FSTmoy=mean(FST,na.rm=T) ; FSTmoy_corr=mean(NUMERATOR)/mean(DENOMINATOR)
  
  list(FST=FST,FSTmoy=FSTmoy,FSTmoy_corr= FSTmoy_corr,num=NUMERATOR,den=DENOMINATOR)
}

###
###
### --- GENOME SCAN ---
###
###

#' BI-COLOUR CHROMOSOMES COLOR FOR MANHATTAN PLOTS
#' 
#' From a list of chromosomes values, it creates a vector of 
#' alternating colors for each chromosome
#' @param list A list of chromosome names
#' @param vector A list of two elements of colors
#' @return A vector of colores of the same size as the list of chromosomes
#' @export
#' 
colorChromosomes <- function(x, colors=c("black", "grey"))
{
  chrm.splitted <- strsplit(x, split = "_")
  chrm.splitted <- as.numeric(do.call(rbind.data.frame, chrm.splitted)[,2])
  
  chrm.splitted.counts <- table(chrm.splitted)
  
  chrm.colors=NULL
  for (i in seq_along(chrm.splitted.counts))
  {
    if((i %% 2) == 0) 
    {
      #print(paste(i,"is Even"))
      c = rep(colors[1],chrm.splitted.counts[i])
    } else {
      #print(paste(i,"is Odd"))
      c = rep(colors[2],chrm.splitted.counts[i])
    }
    chrm.colors <- c(chrm.colors, c)
  }
  
  return(chrm.colors)
}

#' CALCULATE ALLELE FREQUENCY CORRELATION WITHIN A WINDOW
#'
#' Calculate allele frequencies correlation between two populations
#' of SNPs within a window
#' @param dataFrame A dataframe with 'CHROM', 'POS', 'AF_1', 'AF_2' as mandatory fields.
#' @param step The sliding window increment step size.
#' @param windowSize The size of each window.
#' @param method One of "pearson" (default), "kendall", or "spearman": can be abbreviated.

alleleRefFreqCorSlidingWindow <- function(dataFrame, step=1000, windowSize=10000, method="pearson")
{
  if (max(dataFrame[,"POS"]) > windowSize)
  {
    w_cor = NULL
    numOfChunks = (max(dataFrame[,"POS"])-windowSize)/step
    stepStarts  = seq(from=1, to=round(numOfChunks*step), by=step)
    for (i in 1:length(stepStarts))
    {
      idx = (dataFrame[,"POS"] >= stepStarts[i]) & (dataFrame[,"POS"] < stepStarts[i] + windowSize)
      window_cor = cor(dataFrame[idx,"AF_1"], dataFrame[idx,"AF_2"], method=method)
      if (!is.nan(window_cor))
      {
        w_cor[i] = window_cor
      } else {
        w_cor[i] = NA
      }
      
    }
    
    idx = (dataFrame[,"POS"]>=(round(numOfChunks*step))) & (dataFrame[,"POS"] <=round(numOfChunks*step)+windowSize)
    window_cor = cor(dataFrame[idx,"AF_1"], dataFrame[idx,"AF_2"], method=method)
    w_cor = c(w_cor,window_cor)
    stepStarts = c(stepStarts, max(stepStarts)+step)
    output= data.frame(CHROM=dataFrame[1,"CHROM"], Step=stepStarts, windowCor=w_cor)
    
  } else {
    stepStarts = 1
    w_cor = cor(dataFrame[idx,"AF_1"], dataFrame[idx,"AF_2"], method=method)
    output= data.frame(CHROM=dataFrame[1,"CHROM"], Step=stepStarts,  windowCor=w_cor)
  }
  
  return(output)
}

###
###
### --- OTHER FUNCTIONS ---
###
###

#' CALCULATE MULTIPLE CORRELATION STATISTICS BETWEEN ELEMENTS IN A LIST OF MATRIX
#'
#' Compute correlation, mean squared error, R squared and bias
#' between vectors of observed and predicted values organized
#' as two distinct matrix in a list of matrix.
#' @param list The list with two matrices (e.g. observed and estimated AF of each population)
#' @return a data.frame ncol_matrix * 4 (summary statistics)
#' @export
#'
calculate.multiple.correlations <- function(list)
{
  m <- do.call(cbind, list)
  cor = NULL
  mes = NULL
  rsquared = NULL
  bias = NULL
  for (i in 1:(dim(m)[2]/2))
  {
    t = cbind(m[,i], m[, (i + (dim(m)[2]/2))])
    t = t[complete.cases(t), ]
    cor <- c(cor, round(cor(t[,1], t[,2]), 3))
    mes <- c(mes, round(mean((t[,2] - t[,1])^2, na.rm = TRUE), 3))
    rsquared <- c(rsquared, round(1 - (mean((t[,2] - t[,1])^2, na.rm = TRUE)/var(t[,1], na.rm = TRUE)),3))
    bias <- c(bias, round(mean(t[,2] - t[,1], na.rm = TRUE),3))
  }
  
  res <- data.frame(cor=cor, mes=mes, rsquared=rsquared, bias=bias)
  return(res)
}

#' EXPORT SCAFFOLD_POSITION WITH pFST > THRESHOLD
#' @param x data.frame with pairwise FSTs
#' @return void
#' @export
export.pfst <- function(x, threshold=0.25)
{
  x_names <- names(x[, 5:length(colnames(x))])
  for (i in seq_along(x_names))
  {
    dt <- x[, c(1,2,(4+i))]
    s  <- dt[which(dt[,3] > threshold), ]
    l  <- paste0(s[,1], '_', s[,2])
    write.table(l, file = paste0(x_names[i],".txt"), 
                col.names = F, row.names = F, quote = F)
    
  }
}

#' FUNCTION TO CALCULATE FST P-VALUES
#' @param x data.frame with pairwise FSTs
#' @param num1 number of chromossomes in population 1
#' @param num2 number of chromossomes in population 2
#' @param method method to control the type I error
#' @return list
#' @export
calculte.pvalue <- function(x, num1, num2, method="BH")
{
  obs_chisq    <- 2*((num1 + num2)/2)*x #Chisq distribution
  chisq_pvalue <- 1-pchisq(obs_chisq, 1)
  chisq_padj   <- p.adjust(chisq_pvalue, method=method)
  minuslog10   <- -log10(chisq_padj)
  
  return(list(chisq=obs_chisq, pvalue=chisq_pvalue, padj=chisq_padj, logpvalue=minuslog10))
}

#' PLOT OMEGA AFTER SVD
#' 
#' Compute single value decomposition on the omega matrix and return a PC plot.
#' Modified from baypass_utils.R.
#' @param omega The omega matrix returned by BAYPASS core model.
#' @param PC a vector of two elements indicating which PCs to plot.
#' @param pop.names a vector with the pool names.
#' @param main a string with the plot title.
#' @param col a vector with the same size as number of pools with color definition.
#' @return a plot
#' @return a list
#' @export
plot.omega.mod <- function(omega, PC=c(1,2), pop.names=paste0("Pop",1:nrow(omega)), 
                       main=expression("SVD of "*Omega), col=rainbow(nrow(omega)), pch=16,pos=2)
{
  om.svd=svd(omega)
  eig=om.svd$d
  pcent.var=100*eig/sum(eig)
  plot(om.svd$u[,PC],main=main,pch=pch,col=col,
       xlab=paste0("PC",PC[1]," (",round(pcent.var[PC[1]],2),"%)"),
       ylab=paste0("PC",PC[2]," (",round(pcent.var[PC[2]],2),"%)"),
       xlim=c((min(om.svd$u[,PC[1]])-0.1 ), (max(om.svd$u[,PC[1]])+0.1)),
       ylim=c((min(om.svd$u[,PC[2]])-0.1 ), (max(om.svd$u[,PC[2]])+0.1)),
  )
  text(om.svd$u[,PC[1]],om.svd$u[,PC[2]],pop.names,col=col,pos=pos)
  list(PC=om.svd$u,eig=eig,pcent.var=pcent.var)
}

#' MAKE SLURM SCRIPT TO RUN PODS TO CALIBRATE XtX STATISTICS
#' @param n_nodes
#' @param n_cores
#' @param n_hours
#' @param memory
#' @param account
#' @param g_file
#' @param poolsize_file
#' @param outprefix
#' @param nthreads
#' @param npilot
#' @param pilotlength
#' @param burnin
#' @param logfile
#' @param outdir
#' @param wd
#' @return vector containing the run specifications

make_baypassrun_slurm_pods <- function(n_nodes=1, n_cores=5, n_hours=3, memory=10, account='PAS1715',
                                       g_file=g_file, poolsize_file=poolsize_file, outprefix=outprefix,
                                       nthreads=5, npilot=25, pilotlength=500, burnin=2500, logfile=logfile,
                                       outdir=outdir, wd=wd)
{
  
  parm = as.data.frame(cbind(n_nodes=n_nodes, n_cores=n_cores, n_hours=n_hours, account=account,
                             nthreads=nthreads, npilot=npilot, pilotlength=pilotlength, burnin=burnin, logfile=logfile))
  
  run_baypass_sh = 'run_baypass_pods.sh'
  
  sink(file=run_baypass_sh, type='output');
  
  # SLURM HEADER
  cat(
    paste0(
      '#!/usr/bin/env bash', '\n',
      '#SBATCH -J run_baypass_pods', '\n',
      '#SBATCH -c', ' ', n_cores, '\n',
      '#SBATCH -N', ' ', n_nodes, '\n',
      '#SBATCH -t', ' ', n_hours, ':00:00', '\n', 
      '#SBATCH --mem', ' ', memory, 'G', '\n',
      '#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/run_baypass.%A_%a.out # Standard output', '\n',
      '#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/run_baypass.%A_%a.err # Standard error', '\n',
      '#SBATCH --account', ' ', account,
      '\n',
      'module load baypass', 
      '\n',
      'wd="/fs/scratch/PAS1715/aphidpool"', 
      '\n',
      'echo "running baypass"', 
      '\n',
      'baypass -gfile', ' ', g_file, ' ', '-poolsizefile', ' ', poolsize_file, ' ', '-outprefix', ' ', outprefix, ' ',
      '-nthreads', ' ', nthreads, ' ', '-npilot', ' ', npilot, ' ', '-pilotlength', ' ', pilotlength, ' ', '-burnin', ' ', burnin, ' ',
      '>', ' ', logfile,
      '\n',
      'echo "moving files to baypass output folder"',
      '\n',
      'mv', ' ', outprefix,'.log', ' ', outdir, ' ', '\n',
      'mv', ' ', outprefix,'_*', ' ', outdir, ' ',
      '\n',
      'echo "done"'
    )
  )
  sink();
  
  #system(paste0('sinteractive', ' -N ', n_nodes, ' -c ', n_cores, ' -t ', n_hours, ':10:00', ' -J ', 'run_baypass_pods', ' -A ', account))
  #system(paste0('wd=', wd))
  #system(paste0('sbatch ${wd}/', run_baypass_sh))
  return(parm)
  
}


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
