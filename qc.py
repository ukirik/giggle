from Bio import SeqIO
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import sys, csv, argparse, re, math

""" Command-line parser """
parser = argparse.ArgumentParser()
parser.add_argument("gene", help="Gene to extract")
parser.add_argument("imgt", help="IMGT file")
parser.add_argument("fastq", help="FASTQ file")
parser.add_argument("-o", "--outfile", help="name of the output file")
parser.add_argument("-t", "--threshold", type=float, default=99.4, help="Threshold value for extracting sequences")
parser.add_argument("-m", "--motif", help="The motif to match sequences to")
parser.add_argument("-v", "--verbose", action="store_true", help="increase verbosity")
args = parser.parse_args()

""" CONSTANTS """
SEQID = 1
VGENE = 3
VIDEN = 5
SEQIND = 28


def errstr(s):
    return "\x1b[1;31;40m{}\x1b[0m".format(s)


def emphstr(s):
    return "\x1b[6;31;47m{}\x1b[0m".format(s)


def getMotifOptions(motif):
    match = re.search("\[(.*)\]", motif, re.IGNORECASE)
    chars = set(match.group(1))

    options = {}
    for c in chars:
        key = re.sub("\[.*\]", c, motif, flags=re.I)
        options[key] = c
    return options


def getMotifRep(motif):
    match = re.search("\[(.*)\]", motif, re.IGNORECASE)
    chars = match.group(0)
    return motif.replace(chars, "X")


def getEmphMotif(s, size=1):
    l = len(s)
    start = int(l / 2)
    stop = start + size
    return s[:start] + emphstr(s[start:stop]) + s[stop:]


def processIMGT():
    """ Processes IMGT file (1_Summary) and returns
    a set of sequence IDs """

    with open(args.imgt) as f:
        seqids = set()
        reader = csv.reader(f, delimiter="\t")
        regex = re.compile(r"_.*")

        for row in reader:
            # filters
            if len(row) < SEQIND:
                continue

            """ Ugly hack to avoid fishing out wrong genes
            e.g. IGHV1-24 instead of IGHV1-2 """
            if args.gene + "*" not in row[VGENE]:
                continue

            if args.threshold < float(row[VIDEN]):
                continue

            if args.motif is not None:
                if motif_regex.search(row[SEQIND]) is None:
                    continue

            seq_id = regex.sub("", row[SEQID])
            seqids.add(seq_id)

        return seqids


def fastqoutput(seqs):
    """ Prints out the sequences to a FASTQ file """

    if args.outfile is not None:
        nseqs = SeqIO.write(seqs, args.outfile, "fastq")
    else:
        nseqs = SeqIO.write(seqs, "qc_seqeunces.fastq", "fastq")

    print("{0} sequences written to file".format(nseqs))


def extractMotifs(record):
    """ Extracts a subsequences from a given set of sequences """

    matches = motif_regex.finditer(str(record.seq))
    matches = list(matches)
    if len(matches) == 1:
        m = matches[0]
        start = m.start()
        stop = m.end()
        return record[start:stop]
    else:
        print(errstr("Unexpected number of matches for motif in sequence"), file=sys.stderr)
        for m in matches:
            print(errstr("{:3d}-{:3d}: {}".format(m.start(), m.end(), m.group(0))), file=sys.stderr)
        print(errstr("Skipping sequence {}".format(record.id)), file=sys.stderr)
        return


def main():

    ids = processIMGT()
    if len(ids) == 0:
        print(errstr("No ids gathered"), file=sys.stderr)
        sys.exit(0)
    else:
        print("{0} sequence IDs extracted".format(len(ids)))

    records = SeqIO.parse(args.fastq, "fastq")
    filterd = [rec for rec in records if rec.id.split('|')[0] in ids]
    if len(filterd) == 0:
        print(errstr("Nothing left after filtering"), file=sys.stderr)
        sys.exit(0)
    else:
        print("{0} sequences are filtered".format(len(filterd)))

    # fastqoutput(filterd)
    trimd = [extractMotifs(rec) for rec in filterd if rec is not None]
    data = defaultdict(list)
    empty_recs = 0
    for rec in trimd:
        if rec is None:
            empty_recs += 1
            continue
        else:
            data[str(rec.seq)].append(rec.letter_annotations["phred_quality"])

    if empty_recs > 0:
        print(errstr("{} empty records found".format(empty_recs)), file=sys.stderr)

    keys = getMotifOptions(args.motif)
    for k in sorted(data.keys()):
        df = pd.DataFrame(data[k])
        df.columns = list(getMotifRep(args.motif))

        if args.verbose:
            print("Motif: {}, size: {}".format(getEmphMotif(k), df.shape))
            print(df.head())

        lab = "{} (n = {})".format(keys[k], len(df.index))
        err = df.std() / math.sqrt(len(df.index))
        ax = df.mean().plot(label=lab, fmt='o-', grid='on', yerr=err)

        if args.verbose:
            plotdata = pd.concat([df.mean(), err], axis=1)
            plotdata.columns = ["Mean Q", "Std Err"]
            print(plotdata)

    caption = "Quality scores for sequence fragment\n{}".format(args.motif)
    ax.set(title=caption, xlabel="Position", ylabel="PHRED Quality")

    ax.axis(ymin=15, ymax=40)
    plt.legend(loc=0)

    if args.outfile is not None:
        outfile = args.outfile + ".pdf"
    else:
        outfile = args.gene + ".pdf"

    plt.savefig(outfile, bbox_inches='tight', format='pdf')
    # plt.show()

if __name__ == "__main__":
    matplotlib.style.use("ggplot")
    if args.motif is not None:
        args.motif = args.motif.upper()
        motif_regex = re.compile(args.motif, re.IGNORECASE)

    main()
