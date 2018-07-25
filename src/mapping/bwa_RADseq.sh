# This script runs the mapping for all samples using the BWA aligner.
# It takes 3 arguments:
# -mm : the number of mismatches allowed in bwa-aln
# -ref : the path to the indexed reference genome
# -reads : the path to the input reads to map to the genome
# Cyril Matthey-Doret
# 11.10.2017

# LSF utilities
source src/misc/jobs_manager.sh
# Default number of mismatches is 4
MM=4
threads=2
prefix=BWA

while [[ "$#" > 1 ]]; do case $1 in
    # Number of mismatches allowed in BWA-aln
    --mm) MM="$2";;
    # Path to indexed ref genome
    --ref) index="$2";;
    # Path to input fastq files to map
    --reads) data_dir="$2";;
    # Path to output alignment files
    --out) out_dir="$2";;
    # Path to popmap file
    --popmap) pom="$2";;
    # optional: run locally
    --local) local=yes;;
    *) break;;
  esac; shift; shift
done


# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi


## create output directory for bam files:
mkdir -p $out_dir/bam

for sample in $(cut -f1 "$pom") #this is the list of sample names
do
  # Skip iteration for this sample if fastq file does not exist
  if [ ! -f $data_dir/$sample.fq* ]
  then
    continue
  fi
  # Limit number of queued mapping jobs to 100 at a time
  if [ -z ${local+x} ];then bmonitor bam 100; fi
  echo "processing sample $sample";
  # Sending each sample as a separate jobs
  eval $run_fun <<MAPSAMPLE
    #!/bin/bash

    #BSUB -L /bin/bash
    #BSUB -o BWA-%J-OUT.txt
    #BSUB -e BWA-%J-ERROR.txt
    #BSUB -u cmatthey@unil.ch
    #BSUB -J bam
    #BSUB -n 2
    #BSUB -R "span[ptile=2]"
    #BSUB -q priority
    #BSUB -R "rusage[mem=4000]"
    #BSUB -M 4000000

    # Loading softwares
    source src/misc/dependencies.sh

    # align reads
    bwa aln -n $MM -t $threads $index $data_dir/$sample.fq.gz > $out_dir/$sample.sai
    # index alignment file
    bwa samse -n 3 $index $out_dir/$sample.sai $data_dir/$sample.fq.gz > $out_dir/$sample-$prefix.sam

    # perl script removes reads which map more than once
  	perl src/mapping/split_sam.pl -i $out_dir/$sample-$prefix.sam -o $out_dir/$sample-$prefix >> $out_dir/split_summary.log

    # Remove original SAM files
    rm -v $out_dir/$sample-$prefix.sam

    # Convert SAM files to BAM
  	samtools view -@ $threads -bS -o $out_dir/bam/$sample-$prefix-uniq.bam $out_dir/$sample-$prefix-uniq.sam

    # Sort alignments by leftmost coordinate
    samtools sort -@ $threads $out_dir/bam/$sample-$prefix-uniq.bam -o $out_dir/bam/$sample-$prefix-uniq.sorted.bam

    # Index BAM files
    samtools index $out_dir/bam/$sample-$prefix-uniq.sorted.bam

    # Output index statistics
  	samtools idxstats $out_dir/bam/$sample-$prefix-uniq.sorted.bam

    # Remove unsorted bam files
  	rm -v $out_dir/bam/$sample-$prefix-uniq.bam
  	date

    # Compress SAM files
    gzip -v $out_dir/$sample*.sam
MAPSAMPLE
done

# Wait for all mapping jobs to be finished before resuming pipeline
if [ -z ${local+x} ];then bmonitor bam 0; fi
