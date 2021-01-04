import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import math

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "\npython %prog \
--input input.vcf.gz \
--step 10000 \
--window 50000 \
--min-count 2 \
--pool-size 80 \
--min-sites-frac 0.75 \
--BED test.bed \
--sample AT_Mau_14_01 \
--output out"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, """
H E L P:
____________

This script calculates Pool-Seq corrected Population genetic parameters pi, Theta and Tajima's D using the corrections in Kofler et al. 2011 (PLoS One) for a given sample (--sample) in the header of the vcf file input (--input). The script produces output files for each statistic. The path and the prefix of these files need to be provided (--output). Average statistics are calculated in windows, were step-sizes are defined by --step and window-sizes by --window. If you want to calculate non-overlapping windows you need to make sure that window-sizes are of the same size as step-sizes. The number of chromosomes in the pool (--pool-size) needs to be provided to apply the Pool-Seq corrections in Kofler et al. 2011 (PLoS One).

The scripts needs a mandatory input file (--BED) in BED file format, which is needed to adjust the window-size by excluding the number of "improper sites" (sites, which are located within TE's, in the proximity to InDels and which did not fail the quality criteria during SNP calling) from a given window.

Additional optional parameters are "minimum allele count" (--min-count, default=2) where only alleles above this thresholds are considered and "Minimum fraction of covered sites" (--min-sites-frac, default=0.75) where only windows are conssidered whose porportion of "proper sites" (see above) is equal or larger than the threshold.

""")
#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="vcf", help="vcf file with all SNPs")
parser.add_option("--step", dest="stp", help="the step-size, e.g. 100000")
parser.add_option("--window", dest="win", help="the window-size, e.g. 100000")
parser.add_option("--min-sites-frac", dest="mins",
                  help="minimum fraction of proper sites in a given window, default=75%", default=0.75)
parser.add_option("--min-count", dest="minc",
                  help="minimum allele count, default=2", default=2)
parser.add_option("--pool-size", dest="pools",
                  help="Pool-Size of target sample")
parser.add_option("--BED", dest="sites", help="BED file with missing data")
parser.add_option("--sample", dest="sample", help="name of target sample")
parser.add_option("--output", dest="out",
                  help="The path and name used for the output-file(s)")

parser.add_option_group(group)
(options, args) = parser.parse_args()


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


def average(x):
    return sum(x) / float(len(x))


def binom(x, y):
    ''' calculate the binomial coefficient for x over y'''
    import math
    from decimal import Decimal
    if y == x:
        return 1
    elif y == 1:
        return x
    elif y > x:
        return 0
    else:
        a = math.factorial(x)
        b = math.factorial(y)
        c = math.factorial(x - y)
    try:
        div = a / (float(b) * float(c))
    except:
        div = Decimal(a) / (Decimal(b) * Decimal(c))
    return float(div)

################ Theta a la Kofler 2011 #######################


def ThetaCorr(M, n, b):
    ''' as in Kofler et al. 2011; Page 7'''
    T = 0.0
    for m in range(int(b), int(M - b + 1)):
        K = 0.0
        for k in range(1, int(n)):
            BE = binom(M, m)
            T1 = ((k / float(n))**m)
            T2 = ((((n) - k) / float(n))**(M - m))
            K += (1 / float(k)) * BE * T1 * T2
        T += K
    if T == 0:
        return "NA"
    return T

################# Pi a la Kofler et al. 2011 #####################


def pi(x, M):
    ''' calculate pi on a SNP-wise basis. where x is a vector of all allelefreqs and n is the samplesize)
    x is a list of allele frequencies
    n is the total coverage'''
    if M == 0:
        return "NA"
    else:
        corr = (M) / float(M - 1)
        freqsum = sum([y**2 for y in x])
        return (1 - freqsum) * corr


def PiCorr(M, n, b):
    ''' see Kofler et al. 2011; Page 7
    x is a list of allele counts
    n is the poolsize
    b is the minor allele threshold'''
    T = 0.0
    for m in range(int(b), int(M - b + 1)):
        UT = (2 * m * (M - m) / float(M * (M - 1)))
        K = 0.0
        for k in range(1, int(n)):
            BE = binom(M, m)
            T1 = ((k / float(n))**m)
            T2 = ((((n) - k) / float(n))**(M - m))
            K += (1 / float(k)) * BE * T1 * T2
        T += UT * K
    if T == 0:
        return "NA"
    return T

################# Tajima's D a la Kofler et al. 2011 #############
# N= poolsize
# C= coverage list from SNPs


def an(n):
    '''Kofler et al. 2011 Appendix: P.6 Formula 4'''
    toret = 0.0
    for i in range(1, int(n)):
        toret += 1 / float(i)
    return toret


def bn(n):
    '''Kofler et al. 2011 Appendix: P.7 Formula 1'''
    toret = 0.0
    for i in range(1, int(n)):
        toret += 1 / float((i**2))
    return toret


def fstar(n, AN):
    '''Kofler et al. 2011 Appendix: P.6 Formula 1'''
    return (n - 3) / (AN * (n - 1) - n)


def alphas(n):
    '''Kofler et al. 2011 Appendix: P.6 Formula 2'''
    AN = an(n)
    FS = fstar(n, AN)
    t1 = (FS**2) * (AN - (n / (n - 1)))
    st1 = AN * ((4 * (n + 1)) / ((n - 1)**2))
    st2 = 2 * ((n + 3) / (n - 1))
    t2 = FS * (st1 - st2)
    t3 = AN * ((8 * (n + 1)) / (n * ((n - 1)**2)))
    t4 = ((n**2) + n + 60) / (3 * n * (n - 1))
    #print "AS:",t1 + t2 - t3 + t4
    return t1 + t2 - t3 + t4


def betas(n):
    '''Kofler et al. 2011 Appendix: P.6 Formula 3'''
    AN = an(n)
    BN = bn(n)
    FS = fstar(n, AN)

    t1 = (FS**2) * (BN - ((2 * (n - 1)) / ((n - 1)**2)))
    st1 = BN * (8 / (n - 1))
    st2 = AN * (4 / (n * (n - 1)))
    st3 = ((n**3) + 12 * (n**2) - 35 * n + 18) / (n * ((n - 1)**2))
    t2 = FS * (st1 - st2 - st3)
    t3 = BN * (16 / (n * (n - 1)))
    t4 = AN * (8 / ((n**2) * (n - 1)))
    st4 = 2 * (n**4 + 110 * (n**2) - 255 * n + 126)
    st5 = 9 * (n**2) * ((n - 1)**2)
    t5 = st4 / st5
    #print "BS:",t1 + t2 - t3 + t4 + t5
    return t1 + t2 - t3 + t4 + t5


def nbase(np, M):
    pij = pijmatrix(3 * np, np)
    nb = 0.0
    minj = max([M, np])
    for k in range(1, minj + 1):
        nb += k * pij[M][k]
    #print nb
    return nb


def pijmatrix(MC, np):
    from collections import defaultdict as d
    jb = min([MC, np])
    Matrix = d(lambda: d(float))
    #print np
    Matrix[0][0] = 1.0
    for i in range(1, MC + 1):
        for j in range(1, jb + 1):
            t1 = ((1 + np - j) / float(np)) * Matrix[i - 1][j - 1]
            t2 = (j / float(np)) * Matrix[i - 1][j]
            pij = t1 + t2
            Matrix[i][j] = pij

    return Matrix


def div(C, theta, avn):
    import math
    AS = alphas(avn)
    BS = betas(avn)
    #print avn,AS,BS,len(C),theta
    div1 = (AS / C) * theta + BS * (theta**2)
    return math.sqrt(abs(div1))


def D(pi, theta, C, avn):
    '''Kofler et al. 2011 Appendix: P.7 Formula 3'''
    Div = div(C, theta, avn)
    if pi - theta == 0:
        return 0
    else:
        return (pi - theta) / Div

################# sliding windows ##############


def slider(x, window, step, TYP):
    '''return list of items in sliding windows, x must be a list of numbers'''
    steps = {}
    stop = max(x)
    start = 0
    end = window
    get = [y for y in set(x) if y > start and y <= end]
    while(end <= stop):
        if TYP == "count":
            steps[(start + end) / 2.0] = len(get)
        else:
            steps[(start + end) / 2.0] = get
        start = start + step
        end = start + window
        get = [y for y in set(x) if y > start and y <= end]
    if get != []:
        if TYP == "count":
            steps[(start + end) / 2.0] = len(get)
        else:
            steps[(start + end) / 2.0] = get
    return steps

################ read parameters ##########


window = [int(x) for x in options.win.split(",")]
step = [int(x) for x in options.stp.split(",")]
mins = float(options.mins)
mincount = float(options.minc)
sample = options.sample
pools = int(options.pools)

OPl = []
OTl = []
ODl = []

for i in range(len(window)):
    OPl.append(open(options.out + "_" +
                    str(window[i]) + "_" + str(step[i]) + ".pi", "w"))
    OTl.append(open(options.out + "_" +
                    str(window[i]) + "_" + str(step[i]) + ".th", "w"))
    ODl.append(open(options.out + "_" +
                    str(window[i]) + "_" + str(step[i]) + ".D", "w"))

PopGendata = d(lambda: d(lambda: d(float)))
avndict = {}
Pcorrhash = {}
Tcorrhash = {}
MissDat = d(list)
count = 1

# calculate adjusted window-sizes:
for l in load_data(options.sites):
    a = l.rstrip().split()
    for x in range(int(a[1]) + 1, int(a[2]) + 1):
        MissDat[a[0]].append(x)
print("BED file read")
WIN = d(lambda: d(str))
for k, v in MissDat.items():
    print("processing", k
          )
    for W in range(len(window)):
        WIN[k][W] = slider(v, window[W], step[W], "count")

covh = d(lambda: d(lambda: d(int)))

for k, v in WIN.items():
    for W in range(len(window)):
        for win, v1 in v[W].items():
            covh[k][W][win] = window[W] - v1

print("Adjusted Window sizes calculated")

for l in load_data(options.vcf):

    a = l.rstrip().split()

    # skip header
    if l.startswith("##"):
        continue

    # store list of sample names
    if l.startswith("#"):
        header = a[9:]
        continue

    # identify position of sample in list
    if sample in header:
        SP = header.index(sample)
    else:
        print("Sample not in Input file")
        sys.exit()

    # ignore multiallelic sites
    if "," in a[4]:
        continue

    # counter
    if count % 100000 == 0:
        print(count, "SNPs processed")
    count += 1

    # Get allelic information for specific sample
    CODE = a[8].split(":")
    pop = a[9:][SP].split(":")

    POP = {CODE[x]: pop[x] for x in range(len(CODE))}

    if POP["FREQ"] == ".":
        continue

    # calculate pi and theta

    # get frequencies and counts and remove alleles that are not among the two major alleles

    if float(POP["FREQ"]) != 0 and float(POP["FREQ"]) != 1:

        Px = [int(POP["RD"]), int(POP["AD"])]
        px = [1 - float(POP["FREQ"]), float(POP["FREQ"])]

        # get coverage and poolsize
        M = int(POP["DP"])
        n = pools

        # calculate pi
        ID = str(M) + "_" + str(n)
        if ID not in Pcorrhash:
            Pcorrhash[ID] = PiCorr(M, n, mincount)
        if Pcorrhash[ID] == "NA":
            continue
        pix = pi(px, M) / Pcorrhash[ID]

        PopGendata[a[0]][int(a[1])]["P"] = pix
        PopGendata[a[0]][int(a[1])]["M"] = M

        # calculate Theta
        if ID not in Tcorrhash:
            Tcorrhash[ID] = ThetaCorr(M, n, mincount)

        PopGendata[a[0]][int(a[1])]["T"] = 1 / Tcorrhash[ID]

for Chrom, Values in sorted(PopGendata.items()):

    for W in range(len(window)):

        if Chrom not in covh:
            continue

        print(Chrom, "chromosome: calculation started for window-size",
              window[W], "and step-size", step[W])
        # generate a dictionary with bins according to windows and stepsizes containing the corresponding positions
        bins = slider(Values.keys(), window[W], step[W], "values")
        for Bin, Items in sorted(bins.items()):
            if Items == []:
                continue
            print(Bin, "window: calculation started")
            ph = 0.0
            th = 0.0
            mh = []
            # loop through positions in bin
            for i in Items:
                # loop through samples:
                ph += Values[int(i)]["P"]
                th += Values[int(i)]["T"]
                if "M" not in Values[int(i)]:
                    continue
                mh.append(Values[int(i)]["M"])

            # ignore if no covered sites
            if len(mh) == 0:
                TD = "NA"
                Pi = "NA"
                Theta = "NA"

            # test if site count at least min-sites-fraction*windowsize:
            elif covh[Chrom][W][Bin] < window[W] * mins:
                TD = "NA"
                Pi = "NA"
                Theta = "NA"

            else:
                AvCov = int(average(mh))
                ID = str(AvCov) + ":" + str(pools)
                if ID not in avndict:
                    avndict[ID] = nbase(pools, AvCov)

                # calculate average pi and Theta
                Pi = ph / covh[Chrom][W][Bin]
                Theta = th / covh[Chrom][W][Bin]

                # calculate Tajima's D
                TD = D(Pi, Theta, covh[Chrom][W][Bin], avndict[ID])

            # write output
            OPl[W].write(Chrom + "\t" + str(Bin) + "\t"
                         + str(Pi) + "\t" + str(covh[Chrom][W][Bin]) + "\n")
            OTl[W].write(Chrom + "\t" + str(Bin) + "\t"
                         + str(Theta) + "\t" + str(covh[Chrom][W][Bin]) + "\n")
            ODl[W].write(Chrom + "\t" + str(Bin) + "\t"
                         + str(TD) + "\t" + str(covh[Chrom][W][Bin]) + "\n")
