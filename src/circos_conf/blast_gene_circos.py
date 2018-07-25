
# Extracting all BLAST hits between genes (non-self) and parse the results into
# A circos-compatible format.
# Cyril Matthey-Doret
# 24.11.2017
import pandas as pd

#=== SETTING PATH AND VARIABLES ===#
# Tabular output (fmt6) from BLAST
blast = "data/homology/MCScanX/input/MCScanX_in.blast"
# BED file with all genes coordinates
bed = "data/homology/MCScanX/input/MCScanX_genes_conv.bed"
# Significance threshold for E-values
Esig=10**-5
# Output folder
circos_links = 'data/circos/blast.lf.txt'

#=== LOADING FILES ===#
blast_file = pd.read_csv(blast, sep="\t", header=None)
bed_file = pd.read_csv(bed, sep="\t", header=None)
# Selecting relevant columns
bed_file = bed_file.iloc[:,0:4]
# Selecting BLAST results with E-value below threshold, except if self-hit
signif = blast_file.loc[(blast_file[10] < Esig) & (blast_file[0] != blast_file[1])]
# Selecting significant pairs of IDs
out_signif = signif.iloc[:,[0,1]]
# Building circos-compatible table with source + target genes
Q_genes = bed_file.merge(out_signif, left_on = 3, right_on = 0)
QS_genes = bed_file.merge(Q_genes, left_on = 3, right_on = "1_y")
# Excluding irrelevant columns
out_links = QS_genes.iloc[:,[1,2,3,5,6,7]]

# Writing output to file
out_links.to_csv(circos_links, sep = " ", header = False, index = False)
