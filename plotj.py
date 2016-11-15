from itertools import islice
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
parser.add_argument("files", nargs='+', help="tab-separated output of IgDiscover")
parser.add_argument("-o", "--out", help="directory where the files should be output")
parser.add_argument("-t", "--threshold", type=float, default=0.25, help="Percentage of J allele expression for V-association")
parser.add_argument("-s", "--stranslator", help="reference file to translate s-codes")
args = parser.parse_args()

""" Constants: update indices as necessary """
DELIM = "\t"
CIND  = 0
VIND  = 1
DIND  = 2
JIND  = 3
PLTDIR = "jplots"


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

def plotJ(jcounts, pdf):
    d = {}
    for j,count in jcounts.items():
        gene,allele = j.split('*')
        if gene not in d:
            d[gene] = {allele : count}
        else:
            d[gene][allele] = count

    df = pd.DataFrame.from_dict(d)
    df_scaled = df.div(df.sum(axis=0), axis=1)
    print(df_scaled)

    df_scaled.T.plot(kind="barh", stacked=True, legend = False).legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("Proportional allele distribution of IGHJ genes")
    plt.savefig(pdf, bbox_inches='tight', format='pdf')

    # Check if more than 1 IGHJ is dominant
    j_usage = (df_scaled > args.threshold).sum(axis=0)
    dominant_js = list(j_usage[j_usage > 1].index)
    if len(dominant_js) > 0 :
        print("Multiple alleles used for {0}".format(dominant_js))

    
    return dominant_js

def getVJ(j, reader, pdf):
    vj = {}
    jcounter = 0

    for row in reader:
        vcounts = defaultdict(int)
        if args.stranslator is not None:
            s = STranslator(args.stranslator)

        if j in row[JIND]:
            jcounter += 1
            vass = row[VIND] if s is None else s.get(row[VIND])
            jass = row[JIND]
            
            if jass in vj:
                vj[jass][vass] += 1
            else:
                vj[jass] = vcounts
                vj[jass][vass] += 1


    df = pd.DataFrame.from_dict(vj, orient="index")
    df_sorted = pd.DataFrame(df, columns = sorted(df.columns, reverse=True, key = lambda s : sortkey(s)))
    
    # get columns that are below the threshold
    df_filtered = df_sorted[df_sorted.sum(axis=1) > jcounter * args.threshold]

    df_filtered.T.plot.barh(figsize=(8.27, 11.69))
    plt.title("IGHV counts per {0} allele".format(j))
    plt.savefig(pdf, bbox_inches='tight', format='pdf')

def sortkey(genestr):
    noallele = genestr.split('*')[0]
    gfamily, gorder, *crap = noallele.split('-')
    LAST = 99
    if gorder.endswith('D'):
        return (LAST, gfamily)
    else:
        return (int(gorder), gfamily)


def main():

    global VIND, JIND

    print(args.files)

    if len(args.files) > 1:
        check_ncol(args.files)
    
    col_headers(args.files[0]) 
    colinfo = input("By default columns {0} and {1} will be used for V and J genes respectively.\n" \
            "Either confirm or input the desired columns numbers separated by a comma:  'x,y':   ".format(VIND,JIND))

    if colinfo:
        newv, newj = map(int,colinfo.split(','))
        VIND = newv
        JIND = newj

    # Plot style
    matplotlib.style.use("ggplot")
    
    for i in range(len(args.files)):
        f = args.files[i]
        print("Processing file: {0}".format(f))

        jcounts = defaultdict(int)
        if f.endswith(".gz"):
            fh = gzip.open(f, 'rt')
        else:
            fh = open(f, 'r')

        reader = csv.reader(fh, delimiter=DELIM)
        header = next(reader)
        for row in reader:
            jcounts[row[JIND]] += 1
         
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
            dom_js = plotJ(jcounts,pdf)

            if len(dom_js) > 0:
                for j in dom_js:
                    fh.seek(0,0) # reset reader
                    next(reader) # skip headers
                    getVJ(j, reader, pdf)


if __name__ == "__main__":
    main()
