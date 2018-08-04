#!/bin/bash
# This script takes a sorted bam file containing reference-aligned RNA-seq reads
# And the the reference genome in FASTA format as input. It assembles the
# transcriptomes using cufflinks.
# Cyril Matthey-Doret
# 11.11.2017


# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -a bam_file -r ref -o out_folder [-l] [-h]
   -a   sorted bam file containing reference-aligned RNA-seq reads
   -r   reference genome in fasta format
   -o   output folder
   -l   local run. If specified, will not use LSF bsub command
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":a:r:o:lh" opt; do
   case $opt in
   a )  BAM=${OPTARG} ;;
   r )  REF=${OPTARG} ;;
   o )  OUT=${OPTARG};;
   l )  local=yes;;
   h )  usage ;;
   \?)  usage ;;
   esac
done

# Reset output folder if it exists
rm -rf $OUT
mkdir -p $OUT

# If local run, run directly, otherwise submit using bsub
if [ -z ${local+x} ];
then
  run_fun="bsub -K"
else
  run_fun="bash"
fi

eval "$run_fun" <<RNA
#!/bin/bash
#BSUB -J RNA_assembly
#BSUB -q normal
#BSUB -n 36
#BSUB -M 64000000
#BSUB -R "span[ptile=36]"
#BSUB -R "rusage[mem=64000]"

# Loading softwares
source src/misc/dependencies.sh

# Index BAM file
samtools index -b $BAM

# Generate bigwig file with coverage along genome
bamCoverage -b $BAM -o $OUT/coverage.bw

# Assemble transcriptome using cufflinks
cufflinks -p 32 -b $REF -u \
          --library-type ff-firststrand $BAM \
          -o "$OUT"

RNA
