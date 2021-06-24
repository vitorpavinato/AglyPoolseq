#####
### AUXILIARY FUNCTIONS
####

## COMPUTE MAXIMUM LIKELIHOOD IMPUTED SAMPLE ALLELE COUNT
#' Compute the maximum likelihood reference allele count
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