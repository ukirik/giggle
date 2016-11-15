from itertools import islice, accumulate
from operator import itemgetter
from collections import defaultdict
from stranslator import STranslator
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import os, argparse, gzip, csv, pprint

""" Command line parsers """
parser = argparse.ArgumentParser()
parser.add_argument("folders", nargs='+', help="output of IgDiscover runs")
parser.add_argument("-o", "--out", help="directory where the files should be output")
parser.add_argument("-g", "--germline-threshold", type=float, default=0.875, help="Threshold value for considering an allele")
parser.add_argument("-d", "--dcoverage", type=float, default=35, help="Threshold value for D-coverage")
parser.add_argument("-s", "--stranslator", help="reference file to translate s-codes")
args = parser.parse_args()

""" Constants: update indices as necessary """
DELIM = "\t"
CIND  = 0
VIND  = 1
DIND  = 2
JIND  = 3
VERR  = 20
DCOV  = 8
PLTDIR = "haplod"
VFASTA = os.path.join("final","database","human_V.fasta")   
FILTAB = os.path.join("final","filtered.tab.gz")

def check_ncol(files):

    print("Checking %d files for consistency in number of columns..." % len(files))
    ncols = []
    for infile in files:
        with gzip.open(infile,'rt') as f:
            reader = csv.reader(f, delimiter=DELIM)
            row = next(reader)
            ncols.append(len(row))
            if len(ncols) > 1 and ncols[-1] != ncols[-2]:
                sys.exit("Number of columns don't match!: " + str(ncols))

    print("Consistent number of columns found {0}".format(ncols))


def col_headers(infile):
    with gzip.open(infile,'rt') as f:
        reader = csv.reader(f, delimiter=DELIM)
        headers = next(reader)
        print("Column headers for the first file:")
        for i in range(0,len(headers)):
            print("[{0}]: {1}".format(i, headers[i])) 

def getVcandidates(fasta):
    all_vs = []
    het_vs = []
    print("reading fasta file: {0}".format(fasta))
    with open(fasta, 'r') as f:
        counter = 0
        for line in f:
            if line.startswith('>'):
                counter += 1
                genecall, *rest = line.split()
                genecall = genecall.replace('>','')
                gene, allele = genecall.split('*')
                if gene in all_vs:
                    het_vs.append(gene)
                all_vs.append(gene)
                
    print("{0} entries read".format(counter))
    return all_vs, het_vs

def sortkey(genestr):
    noallele = genestr.split('*')[0]
    gfamily, gorder, *crap = noallele.split('-')
    LAST = 99
    if gorder.endswith('D'):
        return (LAST, gfamily)
    else:
        return (int(gorder), gfamily)

def hetVD(tabfile, hetv):
    data = {}
    haplo_vd = defaultdict(lambda: defaultdict(int))

    if args.stranslator is not None:
        s = STranslator(args.stranslator)

    with gzip.open(tabfile, 'rt') as fh:
        reader = csv.reader(fh, delimiter =DELIM)
        header = next(reader)
        for row in reader:
            if not row[DIND] or 'OR' in row[DIND]:  # filter out missing or non-functional D-assignments
                continue
            if float(row[DCOV]) < args.dcoverage:   # filter out D-assignments with little coverage
                continue
            if float(row[VERR]) > 0:                # filter out ambigious V-assignments
                continue
            
            vass = row[VIND] if args.stranslator is None else s.get(row[VIND])
            vgene = vass.split('*')[0]
            dgene = row[DIND].split('*')[0]

            # all V-candidate D pairings
            if vgene in hetv:
                haplo_vd[vass][dgene] += 1

        # filter out uncertain V-alleles 
        for candidate in hetv:
            filterAlleles(haplo_vd, candidate)

        for key,val in haplo_vd.items():
            gene = key.split('*')[0]
            if gene not in data:
                data[gene] = {key: val}
            else:
                data[gene][key] = val

        return haplo_vd, data           

def filterAlleles(haplo_vd, gene):
    temp = {k: sum(v.values()) for k, v in haplo_vd.items() if gene in k}
    sorted_temp = sorted(temp.items(), key = itemgetter(1), reverse=True)

    values = [v for k,v in sorted_temp]
    cumsum = list(accumulate(values))
    scaled = [v / sum(values) for v in cumsum]
    alleles = []
    
    for i in range(0,len(scaled)):
        if scaled[i] > args.germline_threshold:
            j = 0 if i == 0 else i+1    # if gene is not heterozygous remove remaining single allele as well 
            miscalls = sorted_temp[j:]
            for call,count in miscalls:
                del haplo_vd[call]
            return

def plotGene(gene, haplodata, pdf):
    df = pd.DataFrame.from_dict(haplodata, orient="index")
    df_scaled = df.div(df.sum(axis=1), axis=0)
    df_sorted = pd.DataFrame(df_scaled, columns = sorted(df_scaled.columns, reverse=True, key = lambda s : sortkey(s)))
    
    ax = df_sorted.T.plot.barh(figsize=(8.27, 11.69))
    ax.set(xlabel="Normalized frequency")
    plt.title("IGHD combination counts per {0} allele".format(gene))
    plt.savefig(pdf, bbox_inches='tight', format='pdf')


def main():

    print(args.folders)

    # Plot style
    matplotlib.style.use("ggplot")
    
    for i in range(len(args.folders)):
        folder = args.folders[i]

        # Step1: get candidates for heterozygous V-genes
        fastaf = os.path.join(folder, VFASTA)
        allv, hetv = getVcandidates(fastaf)
        print("potentially heterozygous V genes: {0}".format(hetv))
        input("Paused... Press <Enter> to continue... ")

        tabfile = os.path.join(folder, FILTAB)

        # Step3: get VD combination counts
        haplo_vd, alldata = hetVD(tabfile, hetv)
        df = pd.DataFrame.from_dict(haplo_vd, orient="index")
        df_scaled = df.div(df.sum(axis=1), axis=0)
        df_sorted = pd.DataFrame(df_scaled, columns = sorted(df_scaled.columns, reverse=True, key = lambda s : sortkey(s)))
        #print(df_sorted.T)
        #input("Paused... Press <Enter> to conitnue...")
        
        # Check and create folders if necessary
        outdir = PLTDIR
        if args.out is not None:
            if not os.path.exists(args.out):
                os.mkdir(args.out)
            outdir = os.path.join(args.out, PLTDIR)
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        filename = input("Please enter a name for results of this analysis: ")
        while not filename.strip():
            filename = input("Please enter a VALID name for the resultant file: ")


        outfile = os.path.join(outdir, filename+".pdf")
        with PdfPages(outfile) as pdf:
            for gene in sorted(alldata.keys()):
                plotGene(gene, alldata[gene] , pdf)


if __name__ == "__main__":
    main()
