from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import sys, os, argparse, csv, pprint, yaml

""" Command line parsers """
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-c", "--config", help="define a YAML config file for multiple runs")
group.add_argument("-f", "--files", nargs='+', help="output of IMGT/High-VQUEST runs")
parser.add_argument("-g", "--gene", required=True, help="IGHV gene to be analyzed")
parser.add_argument("-p", "--pos", required=True, nargs='+', type=int, help="the particular position, or range of positions, where the SNV occurs")
parser.add_argument("-o", "--out", help="name of the output file")
parser.add_argument("--nofilter", action="store_true", default=False, help="Flag to indicate whether to filter out unproductive sequences (default: %(default)s)")
args = parser.parse_args()

""" Validate pos """
if len(args.pos) > 2:
	sys.exit("Positions argument needs to be a single integer or two integers denoting a range!")

""" Constants: update indices as necessary """
DELIM = "\t"
FUNC = 2
VIND = 3
SIND = 6

def check_cols(files):
	print("Checking %d files for consistency in number of columns..." % len(files))
	col_headers = {}
	for i in range(0, len(files)):
		with open(files[i], 'r') as f:
			reader = csv.reader(f, delimiter=DELIM)
			col_headers[i] = next(reader)
			if len(col_headers) > 1 and col_headers[i-1] != col_headers[i]:
				sys.exit("Column headers don't match!: " + str(ncols))

	for i in range(0, len(col_headers[0])):
		print("[{0}]: {1}".format(i, col_headers[0][i])) 

def get_counts(infile):
	counts = defaultdict(int)
	with open(infile, 'r') as f:
		reader = csv.reader(f, delimiter=DELIM)
		headers = next(reader)
		for row in reader:

			if not args.nofilter and row[FUNC] != "productive": # filter out unproductive (unless filtering is set to false)
				continue

			if args.gene not in row[VIND]: # filter out rows that don't contain the gene of interest
				continue

			if len(args.pos) == 1:
				key = row[SIND][args.pos[0] -1 ] 	# IMGT numbering
			elif len(args.pos) == 2:
				key = row[SIND][args.pos[0] -1 : args.pos[1] -1]
			else:
				print(args.pos, len(args.pos))
				sys.exit("args.pos assertion failed!")

			counts[key] += 1

	return counts

def generatePlot(pdf, files, sample = "N/A"):
	data = {}
	for i in range(0, len(files)):
		# Naming of the series
		# name = "file" + str(i)
		# name = os.path.basename(args.files[i])
		name = os.path.split(os.path.dirname(files[i]))[-1]
		data[name] = get_counts(files[i])

	df = pd.DataFrame.from_dict(data, orient="index")
	df_scaled = df.div(df.sum(axis=1), axis=0)
	df_sorted = pd.DataFrame(df_scaled, index = sorted(df_scaled.index), columns = sorted(df_scaled.columns))
	print(sample)
	print(df_sorted)

	ax = df_sorted.plot.barh(stacked=True, figsize=(8.27, 11.69))
	ax.set(xlabel="Relative frequency", ylabel="Nucleotide(s) at {0}".format(args.pos))
	ax.tick_params(labelsize=10)
	plt.title("Distribution of sequence variation on {0} gene in {1}".format(args.gene, sample))
	plt.legend(loc=9, ncol=5, mode="expand")
	plt.savefig(pdf, bbox_inches='tight', format='pdf')

def main():
	# Plot style
	matplotlib.style.use("ggplot")

	if args.out is None:
		outfile = "SNVplot_" + args.gene + ".pdf"
	else:
		outfile = args.out + ".pdf"

	with PdfPages(outfile) as pdf:

		if args.files is not None:
			check_cols(args.files)
			fig = generatePlot(pdf, args.files)

		elif args.config is not None:
			with open(args.config, 'r') as ymlfile:
				cfg = yaml.safe_load(ymlfile)
				for run in cfg["runs"]:
					generatePlot(pdf, cfg["runs"][run], run)


if __name__ == "__main__":
	main()