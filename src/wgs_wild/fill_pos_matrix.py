# Filling missing position (rows) with NAs in genotype matrix
# Reading input matrix from sys.argv and writing to stdout

from __future__ import print_function
from sys import argv

# Track first row in file
start_line = True
with open(argv[1]) as geno_r:
    for pos in geno_r:
        # Initialize number of samples on first line
        if start_line:
            curr_line = pos.split('\t')
            n_samples = len(curr_line) - 2
            start_line = False
        else:
            prev_line = curr_line[:]
            curr_line = pos.split('\t')
            # Process only if both positions on same chromosome
            if prev_line[0] == curr_line[0]:
                # If gap between previous and current positions, insert lines
                if (int(prev_line[1]) + 1) != int(curr_line[1]):
                    gap_size = int(curr_line[1]) - (int(prev_line[1]) + 1)
                    for i in range(1, gap_size + 1):
                        # Inserting lines with only missing values
                        print('\t'.join([prev_line[0],str(int(prev_line[1]) + i)] + \
                              ['.' for n in range(n_samples)]))
        # Print current line after gap is filled (if needed)
        print('\t'.join(curr_line), end='')
