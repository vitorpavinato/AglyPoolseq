import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import os as OS

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--VCF", dest="VCF", help="Input VCF file")
parser.add_option("--out", dest="out", help="outputfile")
parser.add_option("--private", dest="PRIV", help="Path to csv with private SNPs")
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

SNPhist=d(lambda: d(lambda: d(int)))
out=open(options.out,mode="wt")
PATH=options.PRIV
Priv=d(lambda: d(str))
MAFs=[float(x) for x in options.MAF.split(",")]
Mutl=d(str)
for file in OS.listdir(PATH):
    if not file.endswith(".csv"):
        continue
    File=OS.path.join(PATH, file)
    print(File)
    for l in load_data(File):
        if l.startswith("chr"):
            continue
        a=l.rstrip().split(",")
        Priv[a[0]][int(a[1])]=a[7]

print("Files reading done")

C=1
for l in load_data(options.VCF):
    if C%100000==0:
        print(C,"SNPs read")
    C+=1

    if l.startswith("##"):
        continue
    elif l.startswith("#"):
        header=l.rstrip().split()[9:]
        continue
    a=l.rstrip().split()
    Chr,Pos,N,REF,ALT=a[:5]
    pops=a[9:]
    NW={header[x]:pops[x] for x in range(len(pops))}
    if int(Pos) not in Priv[Chr]:
        continue
    Targetpop=Priv[Chr][int(Pos)]
    for i in range(len(pops)):
        if pops[i].split(":")[0]!="./." and header[i]!=Targetpop:
            AltGT=pops[i].split(":")[0]
            break

    if AltGT=="0/0":
        Mut=REF+">"+ALT
    else:
        Mut=ALT+">"+REF

    if Mut not in Mutl:
        Mutl[Mut]

    for MAF in MAFs:
        if float(NW[Targetpop].split(":")[-1])>=MAF:
            SNPhist[Targetpop][MAF][Mut]+=1

H=""
for k,v in sorted(SNPhist.items()):
    for M,v1 in sorted(v.items()):
        PL=[]
        if H=="":
            out.write("Pop\tMAF\t"+"\t".join(sorted(Mutl.keys()))+"\tCumm\n")
            H=1
        for Mut in sorted(Mutl.keys()):
            if Mut not in v1:
                PL.append(0)
            else:
                PL.append(v1[Mut])
        out.write(str(k)+"\t"+str(M)+"\t"+"\t".join([str(x) for x in PL])+"\t"+str(sum(PL))+"\n")
