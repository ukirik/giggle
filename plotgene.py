from collections import namedtuple, OrderedDict
from itertools import islice
from allele import Allele
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import os, argparse, gzip, csv, pprint

""" Command line parsers """
parser = argparse.ArgumentParser()
parser.add_argument("files", nargs='+', help="tab-separated output of IgDiscover")
parser.add_argument("-o", "--out", help="file to output")
args = parser.parse_args()

""" Constants: update indices as necessary """
DELIM = "\t"
CIND  = 0
VIND  = 1
DIND  = 2
JIND  = 3
VERR  = 20
CDR3  = 28 # This col gives nucleotide seq, use 29 for AA-seq
PLOTDIR = "geneplots"

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

def getAllele(row, d):
    """ Checks if the given row contains an allele that exists already. 
        If so the info for that allele is updated,
        otherwise a new object is created.
    """
    gene = row[VIND].split('*')[0]
    a = Allele(row[VIND])
    if gene not in d:
        d[gene] = [a]
    else:
        alleles = d[gene]
        if a in alleles:
            for allele in alleles:
                if allele == a:
                    allele.addRow(row[CIND], row[DIND], row[JIND], row[CDR3], float(row[VERR]) == 0) 
                    break
        else:
            d[gene].append(a)

def plotgene(genotypes, gene):
    data = {}
    data2= {}
    
    for donor, gt in genotypes.items():
        if gene not in gt:
            continue

        alleles = genotypes[donor][gene]
        datarow = {}
        datarow2= {}
        for a in alleles:
            datarow[a.shortname()] = a.nseqs
            datarow2[a.shortname()]= len(a.cdr3s)
        data[donor] = datarow
        data2[donor]= datarow2

    #donors = ["Donor"+str(i) for i in range(1,len(data)+1)] 
    df = pd.DataFrame.from_dict(data, orient="index")
    df2= pd.DataFrame.from_dict(data2,orient="index")
    #df.index = donors
    #df2.index= donors

    matplotlib.style.use("seaborn-colorblind")
    fig = plt.figure() # Create matplotlib figure

    ax = fig.add_subplot(111) # Create matplotlib axes
    ax2 = ax.twinx() # Create another axes that shares the same x-axis as ax.


    wid = 0.25

    df_scaled = df.div(df.sum(axis=1), axis=0)
    df2_scaled = df2.div(df2.sum(axis=1), axis=0)

    df2_scaled.plot.bar(stacked=True,  ax=ax2, width=wid, position=1,legend=False)
    df_scaled.plot.bar(stacked=True, ax=ax, width=wid, position=0, legend=False)
    fig.autofmt_xdate()
    #ax.set_ylabel("relative sequence count")
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    #ax2.set_ylabel("relative unique CDR3 count")

    #ax.set_ylim(top = 1.5)
    #ax2.set_ylim(top= 1.5)
    plt.title(gene + " expression \n sequence counts (left) and unique CDR3s (right)") 
    
    return plt

def main():

    check_ncol(args.files)
    col_headers(args.files[0])
    genotypes = {}

    for i in range(0,len(args.files)):
        f = args.files[i]
        print("Processing file: {0}".format(f))
        acc = {}
        if f.endswith(".gz"):
            reader = csv.reader(gzip.open(f, 'rt'), delimiter=DELIM)
            header = next(reader)
            for row in reader:
                getAllele(row, acc)
            
            genotypes["Donor"+str(i+1)] = acc
       
    #pprint.pprint(genotypes[0])
    if not os.path.exists(PLOTDIR):
        os.mkdir(PLOTDIR)

    genes = set()
    for donor,gt in genotypes.items():
        for g in gt.keys():
            if g not in genes:
                figure = plotgene(genotypes, g)
                plotfile = os.path.join(PLOTDIR, g+".pdf")
                figure.savefig(plotfile, bbox_inches="tight")
                plt.close()
                genes.add(g)

if __name__ == "__main__":
    main()
