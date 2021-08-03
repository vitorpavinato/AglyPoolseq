#####
### AUXILIARY FUNCTIONS
####

#' COMPUTE MAXIMUM LIKELIHOOD IMPUTED SAMPLE ALLELE COUNT
#'
#' This function obtains the reference allele count with 
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

## COMPUTE HE FROM INPUTED AF
meanHE <- function(x)
{
  l <- sum(!is.na(x))
  return(1 - (sum(((x^2) + ((1-x)^2)), na.rm = T)/l))
}

thetaHE <- function(x)
{
  l <- sum(!is.na(x))
  return((x/(1-x))/l)
}

meanLocusHE <- function(x)
{
  r = 1 - ((x^2) + ((1-x)^2))
  l = sum(!is.na(r))
  return(sum(r, na.rm = T)/(l-1))
}

#' BI-COLOUR CHROMOSOMES COLOR
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

