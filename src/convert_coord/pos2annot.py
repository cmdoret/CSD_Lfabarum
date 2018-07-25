"""
This script takes a genomic position as input (chr,BP), a gff file and a
file with containing GO terms based on the same assembly and returns either \
the annotations within a given region or the top N closest annotations.

Cyril Matthey-Doret
20.10.2017
"""

from __future__ import print_function  # Use python3 print function
import sys  # Will allow to print messages to stderr
import argparse  # Parses command line arguments
import pybedtools as bed  # Python wrapper for bedtools
import pandas as pd
import re

# Parsing command line arguments

parser = argparse.ArgumentParser(description="Takes a genomic position in the \
                                 form 'Chr BP' where Chr is the name of the \
                                 contig and BP is an integer representing the \
                                 basepair position within the contig. Returns \
                                 the annotations in the neighborhood of that \
                                 position. A GFF files with annotations must \
                                 be provided.")
parser.add_argument('--pos', type=str, nargs='?', default=sys.stdin,
                    help='Chromosome name and basepair position of the query. \
                    Should be in the form: Chr,bp where Chr and bp \
                    are a string and integer, respectively. Read from stdin \
                    by default.')
parser.add_argument('gff_file', type=str, help="The path to the GFF file \
                    containing the features.")
parser.add_argument('annot_file', type=str, help="The path to the file containing \
                    the annotations matching feature IDs from the GFF file.")
parser.add_argument('--method', type=str, choices=['range', 'top'],
                    default='range', help='Method used to find annotations, \
                    can be either "range", in which case all annotations \
                    within a basepair range will be picked up, or "top", in \
                    which case the top N nearest annotations will be used. \
                    Default: range. NOTE: only range is available atm.')
parser.add_argument('--range_size', type=int, default=10000, help='In case the \
                    method used it "range", an integer defining the range in \
                    which to look for annotations around the input position.')
parser.add_argument('--top_count', type=int, default=10, help='In case the \
                    method used is "top", the number of closest neighbour \
                    annotations that should be returned.')

args = parser.parse_args()


def eprint(*args, **kwargs):
    # This function prints its arguments to the standard error
    print(*args, file=sys.stderr, **kwargs)


# Handling query position, whether it is from CL arg or stdin
try:
    pos = args.pos.read()
except AttributeError:
    pos = args.pos
# Splitting query into chromosome and basepair
if len(pos.split(',')) == 2:
    Chr = pos.split(',')[0]
    try:
        bp = int(pos.split(',')[1])
    except ValueError:
        eprint("Error: bp must be an integer number representing the \
              genomic position in basepairs within the chromosome.")
elif len(pos.split(',')) == 0:
    eprint("Error: No query provided, please provide a genomic position in the \
           form 'Chr,bp'.")
else:
    eprint("Error: Query position must be in the form 'Chr,bp'.")
    exit()


# Setting up search scope
qstart = max(0, (bp - args.range_size))
qend = bp + args.range_size
# q_region = bed.BedTool('{0} {1} {2}'.format(Chr, qstart, qend),
#                        from_string=True)[0]
eprint("Looking for annotations on chromosome {0}, between \
positions {1} and {2}".format(Chr, qstart, qend))

# Reading GFF and annotations files
gff = bed.BedTool(args.gff_file)
annot = pd.read_csv(args.annot_file, sep='\t', header=None)
annot.columns = ['ID', 'GO', 'term']

# Filter features in query range
gff_filtered = gff.filter(lambda b: (b.chrom == Chr) &
                                    (b.start >= qstart) &
                                    (b.end <= qend))
# Regex matching feature ID in GFF file
feature_pattern = re.compile(r'ID=([^:;-]*)')
# Iterating over GFF features in query region
annot_frame = pd.DataFrame()
for f_id in gff_filtered:
    # Extracting ID
    f_match = re.search(feature_pattern, f_id[8])
    # Retrieving matching annotations
    try:
        f_annot = annot.loc[annot.ID.str.contains(f_match.group(1)), :]
        if not f_annot.empty:
            f_annot.insert(0, 'end', f_id[4])
            f_annot.insert(0, 'start', f_id[3])
            f_annot.insert(0, 'chrom', f_id[0])
            annot_frame = annot_frame.append(f_annot)
    except AttributeError:
        eprint("{0} did not match any GO annotation.".format(f_id))

if not annot_frame.empty:
    annot_frame.drop_duplicates(subset='ID', inplace=True)
    annot_frame.to_csv(sys.stdout, sep='\t', index=False, header=False)
