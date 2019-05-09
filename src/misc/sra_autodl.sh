#!/bin/bash
# Downloads fastq files matching a range of identifiers on SRA
# Requires the fasterq-dump from NCBI
# cmdoret, 20190509

sra_start=$1
sra_end=$2
out_dir=$3

# Extract longest common prefix between identifiers
sra_prefix=$(grep -zPo '(.*).*\n\K\1' <<< "$sra_start"$'\n'"$sra_end")
sra_num_start=${sra_start#????}
sra_num_end=${sra_end#????}
for sra_id in ${sra_prefix}{$sra_num_start..$sra_num_end}; do
        fasterq-dump 
done
