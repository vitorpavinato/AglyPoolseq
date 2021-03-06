import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import gzip

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'simple script that filters by MAF per sample')

#########################################################   CODE   #########################################################################

parser.add_option("--VCF", dest="VCF", help="Input VCF file")
parser.add_option("--MAF", dest="MAF", help="minor allele frequency (MAF) threshold per sample")

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

MAF=float(options.MAF)

for l in load_data(options.VCF):
    if l.startswith("#"):
        print(l.rstrip())
        continue
    a=l.rstrip().split()
    pops=a[9:]
    #print(pops)
    FL=[]
    FR=0
    for i in pops:
        GT,RD,AD,DP,FREQ=i.split(":")
        if GT=="./.":
            FL.append(i)
            continue
        if float(FREQ)<MAF:
            FL.append("./.:.:.:.:.")
            continue
        FL.append(i)
        FR+=1
    if FR==0:
        continue
    print("\t".join(a[:9])+"\t"+"\t".join(FL))
