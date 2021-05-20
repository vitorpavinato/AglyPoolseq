import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import gzip

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Input file")
parser.add_option("--MAF", dest="MAF",
                  help="GLOBAL minimum allele frequency thresholds as comma separated list")
parser.add_option("--Parse", dest="PA", help="parse name", action="store_true")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


if options.PA:
    print(options.IN.split("aphidpool.PoolSeq.PoolSNP.")
          [1].split(".paramTest.ann.vcf.gz"))
    MAF, MAC = options.IN.split("aphidpool.PoolSeq.PoolSNP.")[1].split(
        ".paramTest.ann.vcf.gz")[0].split(".")

stepsize = [float(x) for x in options.MAF.split(",")]
Counts = d(lambda: d(lambda: d(lambda: d(int))))
C = 1
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    if l.startswith("#"):
        a = l.rstrip().split()
        header = a[9:]
        continue
    a = l.rstrip().split()
    if "," in a[4]:
        continue
    INFO = a[7]
    pops = a[9:]
    if C % 100000 == 0:
        print(C, "SNPS processed")
    C += 1
    for i in range(len(pops)):
        if pops[i].startswith("./.") or pops[i].startswith("0/0") or pops[i].startswith("1/1"):
            continue
        for x in stepsize:
            if float(pops[i].split(":")[-1]) < x:
                continue
            if "missense_variant" in INFO:
                Counts[x][i][a[0]]["NS"] += 1
            if "synonymous_variant"in INFO:
                Counts[x][i][a[0]]["SS"] += 1
if options.PA:
    out = gzip.open(options.OUT, "wt")
    out.write("MAF\tMAC\tMAFc\tPOP\tChrom\tNS\tSS\tpNpS\n")
    for MAFc, v in sorted(Counts.items()):
        for Pop, v1 in sorted(v.items()):
            ns, ss = 0, 0
            for chrom, v2 in sorted(v1.items()):
                if v2["SS"] == 0:
                    out.write(MAF[:1] + "." + MAF[1:] + "\t" + MAC + "\t" + str(MAFc) + "\t" +
                              header[Pop] + "\t" + chrom + "\t" + str(v2["NS"]) + "\t" + str(v2["SS"]) + "\tNA\n")
                else:
                    out.write(MAF[:1] + "." + MAF[1:] + "\t" + MAC + "\t" + str(MAFc) + "\t" + header[Pop] + "\t" +
                              chrom + "\t" + str(v2["NS"]) + "\t" + str(v2["SS"]) + "\t" + str(v2["NS"] / v2["SS"]) + "\n")
                ns += v2["NS"]
                ss += v2["SS"]
            if ss == 0:
                out.write(MAF[:1] + "." + MAF[1:] + "\t" + MAC + "\t" + str(MAFc) + "\t" +
                          header[Pop] + "\tgenomewide\t" + str(ns) + "\t" + str(ss) + "\tNA\n")
            else:
                out.write(MAF[:1] + "." + MAF[1:] + "\t" + MAC + "\t" + str(MAFc) + "\t" + header[Pop] +
                          "\tgenomewide\t" + str(ns) + "\t" + str(ss) + "\t" + str(ns / ss) + "\n")
else:
    out = gzip.open(options.OUT, "wt")
    out.write("MAF\tPOP\tChrom\tNS\tSS\tpNpS\n")
    for MAF, v in sorted(Counts.items()):
        for Pop, v1 in sorted(v.items()):
            ns, ss = 0, 0
            for chrom, v2 in sorted(v1.items()):
                if v2["SS"] == 0:
                    out.write(str(MAF) + "\t" + header[Pop] + "\t" + chrom + "\t" + str(
                        v2["NS"]) + "\t" + str(v2["SS"]) + "\tNA\n")
                else:
                    out.write(str(MAF) + "\t" + header[Pop] + "\t" + chrom + "\t" + str(
                        v2["NS"]) + "\t" + str(v2["SS"]) + "\t" + str(v2["NS"] / v2["SS"]) + "\n")
                ns += v2["NS"]
                ss += v2["SS"]
            if ss == 0:
                out.write(str(
                    MAF) + "\t" + header[Pop] + "\tgenomewide\t" + str(ns) + "\t" + str(ss) + "\tNA\n")
            else:
                out.write(str(MAF) + "\t" + header[Pop] + "\tgenomewide\t" + str(
                    ns) + "\t" + str(ss) + "\t" + str(ns / ss) + "\n")
