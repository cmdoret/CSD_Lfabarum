#==================== PARAMETERS ====================#

## Mapping parameters for bwa
# aln: mismatches (MM)
MM=4
## STACKS parameters
# pstacks: minimum coverage (M)
M=3
# cstacks: loci mismatches (LM)
LM=3
# populations: Proportion of individuals with locus (R), Minimum stack depth (D)
R=80
D=5
# Ploidy: exclude haploid males
# Note: 0.90 is a good value to separate ploidy in our dataset (i.e. between the two modes)
# # Change it according to distribution of homozygousity in the dataset
HOM_PLOID=0.90

#======================= PATH =======================#

## Base directories. Avoid changing these.
DAT=./data
LOG=$(DAT)/logs
PROC=$(DAT)/processed

## Scripts
BWA-RAD=src/mapping/bwa_RADseq.sh
P-SRC=src/stacks_pipeline/pstacks_script.sh
C-SRC=src/stacks_pipeline/cstacks_script.sh
S-SRC=src/stacks_pipeline/sstacks_script.sh
GR-SRC=src/stacks_pipeline/group_sstacks.sh
POP-SRC=src/stacks_pipeline/populations_script.sh
RNA-SRC=src/rna_seq/process_rna.sh
MCSX-SRC=src/homology/MCScanX_prep.sh
ASSOC-SRC=src/assoc_mapping/


#======================= DATA =======================#

# Genome-related files
REF=$(DAT)/ref_genome/ordered_genome/merged.fasta
REF-ANN=$(DAT)/ref_genome/ordered_genome/merged.fasta.ann
SIZES=$(DAT)/ref_genome/ordered_genome/tig.sizes.txt
# Mapping data
MAP=$(DAT)/mapped/aln-$(MM)
# STACKS folders
PSTACK=$(DAT)/pstacks/covmin-$(M)
USTACK=$(DAT)/ustacks/cov-$(M)_mm-$(MM)
CSTACK=$(DAT)/cstacks/mm-$(LM)
SSTACK=$(DAT)/sstacks
POP=$(DAT)/populations/grouped_d-$(D)_r-$(R)
# Ploidy-related files
VCFSUM=$(DAT)/ploidy/vcftools/summary_full.txt
THRESH=$(DAT)/ploidy/thresholds/fixed.tsv
BLACK=$(DAT)/ploidy/blacklist.tsv
# Association mapping files
ASSOC=$(DAT)/assoc_mapping
HITS=data/assoc_mapping/case_control/case_control_hits.tsv
# Centromeres files
CENTRO=$(ASSOC)/centro
# RNA-related files (for circos plot)
RNA=$(DAT)/rna_seq/
BAM=$(RNA)/STAR_larvae_aligned.sort.bam.bam
OLD-REF=$(DAT)/ref_genome/canu2_low_corrRE.contigs.fasta
# MCScanX-related files (for circos plot)
MCSX-IN=$(DAT)/homology/MCScanX/input
CORRESP=$(DAT)/annotations/corresp_gff.csv
# Linkage map
LINKMAP=$(DAT)/linkage_map/
# WGS data (used for PI diversity)
WGS=$(DAT)/wgs/
MISC=src/misc
