import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--VCF", dest="VCF", help="Input VCF file")
parser.add_option("--out", dest="out", help="outputfile")
parser.add_option("--MAFs", dest="MAF", help="Minor allele frequency (MAF) thresholds, comma separated")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y

SNPhist=d(lambda: d(int))
MAFs=[float(x) for x in options.MAF.split(",")]
out=open(options.out,mode="wt")

C=1
for l in load_data(options.VCF):
    if C%100000==0:
        print(C,"SNPs read")
    C+=1

    if l.startswith("#"):
        continue

    a=l.rstrip().split()
    pops=a[9:]
    #print(pops)
    FL=[]
    for i in pops:
        GT,RD,AD,DP,FREQ=i.split(":")
        if GT=="./." or max([float(x) for x in FREQ.split(",")])==0:
            continue
        FL.append(max([float(x) for x in FREQ.split(",")]))
    for MAF in MAFs:
        SNPhist[len([x for x in FL if x >=MAF])][MAF]+=1

out.write("POPs\t"+"\t".join([str(x) for x in MAFs])+"\n")
for k,v in sorted(SNPhist.items()):
    PL=[]
    for MAF in MAFs:
        if MAF not in v:
            PL.append(0)
        else:
            PL.append(v[MAF])
    out.write(str(k)+"\t"+"\t".join([str(x) for x in PL])+"\n")
