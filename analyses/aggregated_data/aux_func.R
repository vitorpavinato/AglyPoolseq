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
#' a maximum likelihood imputation appraoch.
#' @param x The the poolseq data(pooldata)
#' @return a matrix of allele frequency counts for each pool in x
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

## COMPUTE GENOME-WIDE MEAN HE FOR EACH SAMPLE
meanHE <- function(x)
{
  l <- sum(!is.na(x))
  return(1 - (sum(((x^2) + ((1-x)^2)), na.rm = T)/l))
}

## COMPUTE GENOME-WIDE MEAN THETA FOR EACH SAMPLE
thetaHE <- function(x)
{
  l <- sum(!is.na(x))
  return((x/(1-x))/l)
}

## COMPUTE HETEROZYGOSITY ACROSS SAMPLES
#'
#' Given the matrix with the ML estimated referece allele frequency
#' for each sample and loci, it computs the locus HE across samples
#' @param x The the ML allele frequency matrix
#' @return a vector with locus HE estimates
#' @export

locusHE <- function(x)
{
  r <- mean(x, na.rm = T)
  a <- 1-r
  h <- 1 - ((r^2) + (a^2))
  return(h)
}

## COMPUTE WITHIN SAMPLE LOCUS GENETIC DIVERSITY - PI OR WITHIN HE
#'
#' It calculates the within population pi (aka HE) for each locus given its imputed
#' reference allele frequencie.
#' @param x The scalar or vector containing the ML imputed reference allele frequencie
#' @param pool_size The 2x the number of individuals in each pool
#' (it corresponds to the number of sampled chromosomes)
#' @return The scalar or vector of within HE
#' @export

withinLocusPi <- function(x, pool_size = 10)
{
    r = x
    a = 1-x
    pi <- (pool_size/(pool_size - 1)) * (r * a * 2)
    return(pi)
}

## COMPUTE LOCUS HETEROZYGOSITY BETWEEN SAMPLES - H1
#'
#' Given the average across samples of the within sample locus genetic diversity
#' and the locus-specific F_ST, it outputs the heterozygosity between samples
#' @param x The scalar or vector containing average across samples of the within
#' sample locus genetic diversity
#' @param locus_fst The scalar or vector of intra-locus F_ST
#' @return The scalar or vector of H1
#' @export

locusH1 <- function(x, f_st = 0.001)
{
    return(x/(1 - (f_st)))
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
