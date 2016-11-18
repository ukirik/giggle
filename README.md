## Description

giggle (or GIgGle) is a collection of Python scripts intended to complement inference tools like TIgGER or IgDiscover. The end goal is to facilitate the interpretation of potential germline variability.

The scripts do have a considerable amount of overlap due to evolving specifications during development. They will be put into a proper modular structure in (hopefully) near future.

## Dependencies
The scripts rely on Pandas, NumPy and matplotlib for visualization. Otherwise it's vanilla Python 3.

## Usage

### haplod.py

Takes in a tab-separated text file(s) as input (by default from IgDiscover) and generates plots based on potential haplotypes, checking associations between which IGHV genes tend to get combined with which IGHD gene. Output figures by default end up in a folder called `haplod` in the current directory, altough it is possible to define a custom folder, rather than the current working directory, by using the `-o flag. 


    usage: haplod.py [-h] [-o OUT] [-g GERMLINE_THRESHOLD] [-d DCOVERAGE]
                     [-s STRANSLATOR]
                     folders [folders ...]

    positional arguments:
      folders               output of IgDiscover runs

    optional arguments:
      -h, --help            show this help message and exit
      -o OUT, --out OUT     directory where the files should be output
      -g GERMLINE_THRESHOLD, --germline-threshold GERMLINE_THRESHOLD
                            Threshold value for considering an allele
      -d DCOVERAGE, --dcoverage DCOVERAGE
                            Threshold value for D-coverage
      -s STRANSLATOR, --stranslator STRANSLATOR
                            reference file to translate s-codes

### plotd.py

Takes in one, or more, tab-separated text file(s) as input (by default from IgDiscover) and generates plots based on usage of IGHD genes with respect the IGHV genes. Output figures by default end up in a folder called `dplots` in the current directory, altough it is possible to define a custom folder, rather than the current working directory, by using the `-o flag. 


    usage: plotd.py [-h] [-o OUT] [-t THRESHOLD] [-d DCOVERAGE] [-s STRANSLATOR]
                    files [files ...]

    positional arguments:
      files                 tab-separated output of IgDiscover

    optional arguments:
      -h, --help            show this help message and exit
      -o OUT, --out OUT     directory where the files should be output
      -t THRESHOLD, --threshold THRESHOLD
                            Percentage of D allele expression for V-association
      -d DCOVERAGE, --dcoverage DCOVERAGE
                            Threshold value for D-coverage
      -s STRANSLATOR, --stranslator STRANSLATOR
                            reference file to translate s-codes
### vdist.py

Takes in one or more IMGT/High-VQUEST output files (specifically the IMGT-gapped-nt-sequences file), looks for a particular position (or a range) in sequences associated with a particular gene, and outputs stacked bar charts visualizing relative distribution of the different nucleotides.

    usage: vdist.py [-h] (-c CONFIG | -f FILES [FILES ...]) -g GENE -p POS
                    [POS ...] [-o OUT] [--nofilter]

    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            define a YAML config file for multiple runs
      -f FILES [FILES ...], --files FILES [FILES ...]
                            output of IMGT/High-VQUEST runs
      -g GENE, --gene GENE  IGHV gene to be analyzed
      -p POS [POS ...], --pos POS [POS ...]
                            the particular position, or range of positions, where
                            the SNV occurs
      -o OUT, --out OUT     name of the output file
      --nofilter            Flag to indicate whether to filter out unproductive
                            sequences (default: False)


### plotgene.py

Takes in one, or more, tab-separated text file as input (by default from IgDiscover) and generates plots based on usage of IGHV genes with respect to sequence counts and number of unique CDR3s. Output figures by default end up in a folder called `geneplots` in the current directory.


### plotj.py

Takes in one, or more, tab-separated text file(s) as input (by default from IgDiscover) and generates plots based on usage of IGHJ genes with respect the IGHV genes. Output figures by default end up in a folder called `jplots` in the current directory, altough it is possible to define a custom folder, rather than the current working directory, by using the `-o flag. 

    usage: plotj.py [-h] [-o OUT] [-t THRESHOLD] [-s STRANSLATOR]
                    files [files ...]

    positional arguments:
      files                 tab-separated output of IgDiscover

    optional arguments:
      -h, --help            show this help message and exit
      -o OUT, --out OUT     directory where the files should be output
      -t THRESHOLD, --threshold THRESHOLD
                            Percentage of J allele expression for V-association
      -s STRANSLATOR, --stranslator STRANSLATOR
                            reference file to translate s-codes

### allele.py

A utility class for abstraction of a IGHV allele, keeping track of number of occurrences, associated -D and -J genes.

### stranslator.py

A utility class for translating `_Sxxxx` codes designated by IgDiscover. Note that this class requires a tabular file which denotes the appropriate name for each `_Sxxxx` pattern, these will be dependent on the original database defined at the time of running IgDiscover. Such a file can be created by submitting the resultant sequence database in the `final` folder of an IgDiscover experiment to IMGT/HighV-QUEST and exporting the results. 



