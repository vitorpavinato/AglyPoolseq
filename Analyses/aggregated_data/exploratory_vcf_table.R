#############################################################
#          VCF parameters: exploratory analysis             #
#                                                           #
#############################################################


## Vitor Pavinato
## vitor.pavinato@supagro.fr
## CFAES

setwd("/fs/scratch/PAS1715/aphidpool")
renv::restore()
## Remove last features
rm(list=ls())
ls()

header <- c("CHROM", "POS", "REF,ALT","ADP","DP","NC", "AC", "AF",
            paste0("DP_", 1:21), paste0("FREQ_", 1:21))

dt <- read.table(file = "vcf/aphidpool.PoolSeq.PoolSNP.001.5.10May2021.vcf.biallele.txt",
                 header = TRUE, col.names = header, sep = "\t", na.strings = ".")

dim(dt)

# Average per-sample depth
hist(dt$ADP, breaks = 100)
min(dt$ADP) #5.18182
max(dt$ADP) #1051.77

(dim(dt)[1]) - dim(dt[dt$ADP < 10, ])[1]
dim(dt[dt$ADP > 500, ])[1]

mean(dt[dt$ADP > 5 & dt$ADP < 10, 6])
plot(table(dt[dt$ADP > 5 & dt$ADP < 10, 6]))
plot(table(dt[dt$ADP > 5, 6]))

mean(dt[dt$ADP > 10 & dt$ADP < 20, 6])
plot(table(dt[dt$ADP > 10 & dt$ADP < 20, 6]))
plot(table(dt[dt$ADP > 10, 6]))

mean(dt[dt$ADP > 20 & dt$ADP < 50, 6])
mean(dt[dt$ADP > 50 & dt$ADP < 100, 6])
mean(dt[dt$ADP > 100, 6])

# DP
hist(dt$DP)
min(dt$DP)
max(dt$DP)
mean(dt$DP)

# AC
hist(dt$AC, breaks = 100)
min(dt$AC)
max(dt$AC)
mean(dt$AC)

mean(dt[dt$AC > 5 & dt$ADP < 10, 6])
length(dt[dt$AC == 10 & dt$ADP < 10,  6])

# AF
hist(dt$AF, breaks = 100)
min(dt$AF)
max(dt$AF)
mean(dt$AF)

plot(dt$AF, dt$AC)
length(dt[dt$AF > 0.01 & dt$AC  == 5, 6])

subdt <- dt[,9:29]

barplot(100*(colSums(is.na(subdt))/265729))
abline(h=30,col="red")

hist(100*(rowSums(is.na(subdt))/21), breaks = 5)

barplot(100*(colSums(is.na(subdt[,-c(5,10,21)]))/265729))


subdt30 <- subdt[100*(rowSums(is.na(subdt))/21) <=30, ]
barplot(100*(colSums(is.na(subdt30))/169721))
abline(h=30,col="red")


barplot(100*(colSums(is.na(dt[dt$ADP > 10, 9:29]))/169363))
abline(h=30,col="red")

hist(100*(dt[dt$NC <= 7 & dt$AC > 10, 6]/21), breaks = 5)

# Using all defined filters

dt.filtered <- dt[dt$NC <= 6 & dt$AC >= 5 & dt$AF >= 0.01, ]
dim(dt.filtered)

barplot(100*(colSums(is.na(dt.filtered[, 9:29]))/151656 ))
abline(h=30,col="red")

# Parameters:
# min cov = 5
# max cov = 0.95
# MAC = 10
# MAF = 0.01
# %MISS = 0.3

# rc (poolfstat) = 2 Make a tri-allelic to a bi-allelic SNP by removing the third allele if AC <= 2


dt.filtered.freq <- dt.filtered[,30:50]

colMeans(dt.filtered.freq, na.rm = T)
colMeans(dt[,30:50], na.rm = T)


pool.names <- c("MN_BIO1_S1",
                "MN_BIO4_S1",
                "ND_BIO1_S1",
                "ND_BIO4_S1",
                "NW_BIO1_S1", # remove
                "NW_BIO1_S2",
                "NW_BIO4_S1",
                "NW_BIO4_S2",
                "NW_BIO4_S3",    
                "PA_BIO1_S1", # remove   
                "PA_BIO4_S1",    
                "PA_BIO4_S2",    
                "WI_BIO1_S1", # remove?    
                "WI_BIO1_S2",    
                "WI_BIO4_S1",    
                "WI_BIO4_S2",    
                "WI_BIO4_S3",    
                "WO_BIO1_S1",    
                "WO_BIO4_S1",    
                "WO_BIO4_S2",    
                "WO_BIO4_S3") # remove
