from itertools import islice
import gzip, argparse, csv, os, sys

""" Command line parsers """
parser = argparse.ArgumentParser()
parser.add_argument("--files", nargs='+', help="tab-separated output of IgDiscover")
parser.add_argument("--gene", nargs='+', required=True, help="IGHV gene to be analyzed")
parser.add_argument("--nofilter", action="store_true", default=False, help="Whether to filter out mutated V-sequences (default: %(default)s)")
parser.add_argument("-s", "--stranslator", help="reference file to translate s-codes")
args = parser.parse_args()

""" Constants: update indices as necessary """
DELIM = "\t"
CIND  = 0
VIND  = 1
DIND  = 2
JIND  = 3
VERR  = 20
CDR3  = 28 	# This col gives nucleotide seq, use 29 for AA
VSEQ  = 30	# This col gives nucleotide seq, use 31 for AA

def check_cols(files):
	print("Checking %d files for consistency in number of columns..." % len(files))
	col_headers = {}
	for i in range(0, len(files)):
		with gzip.open(files[i], 'rt') as f:
			reader = csv.reader(f, delimiter=DELIM)
			col_headers[i] = next(reader)
			if len(col_headers) > 1 and col_headers[i-1] != col_headers[i]:
				sys.exit("Column headers don't match!: " + str(ncols))

	for i in range(0, len(col_headers[0])):
		print("[{0}]: {1}".format(i, col_headers[0][i])) 

def main():

    check_cols(args.files)

    if args.stranslator is not None:
        s = STranslator(args.stranslator)

    for i in range(0, len(files)):
		with gzip.open(files[i], 'rt') as f:
			reader = csv.reader(f, delimiter=DELIM)
			header = next(reader)
			for row in reader:
				


if __name__ == "__main__":
    main()