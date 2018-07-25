#!/bin/bash

# This script runs one pstacks job per sample in a loop on an LSF
# cluster. This speeds up the process by heavily parallelizing tasks.
# Cyril Matthey-Doret
# 10.10.2017

source src/misc/jobs_manager.sh
declare -i ID=1 ## the stacks ID to start from
# Default minimum stack coverage: 3
M=3

# parsing CL arguments
while [[ "$#" > 1 ]]; do case $1 in
    # Mimimum stack coverage for pstacks
    --m) M="$2";;
    # Path to the folder containing input alignment files
    --map) map="$2";;
    # Path to output stacks files
    --out) out_dir="$2";;
    # Path to log folder
    --log) logs="$2";;
    # optional: run locally
    --local) local=yes;;
    *) break;;
  esac; shift; shift
done

# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi

# Cleaning directories
rm -rf $out_dir; mkdir -p $out_dir
rm -rf $logs/pstacks; mkdir -p $logs/pstacks

# For each mapped sample
for i in $(ls $map/bam/*uniq.sorted.bam)
do
    # Do not queue more than 100 jobs at a time
    if [ -z ${local+x} ];then bmonitor PST 100; fi
    echo "Sample= $i, ID=$ID"
    j=$(echo ${i##*/} | cut -f1 -d '.')

    eval $run_fun <<PST
    #!/bin/bash
    #BSUB -L /bin/bash
    #BSUB -o $logs/pstacks/PST_COVMIN${M}_${j}_STDOUT.log
    #BSUB -e $logs/pstacks/PST_COVMIN${M}_${j}_STDERR.log
    #BSUB -J PST${j}
    #BSUB -n 3
    #BSUB -M 2000000
    #BSUB -q priority

    # Loading softwares
    source src/misc/dependencies.sh
    pstacks -f $i -i $ID -o $out_dir -m $M -p 3 -t bam
PST
    ID+=1
done

# Waiting until all Pstacks jobs are finished
if [ -z ${local+x} ];then bmonitor PST 0; fi

# Cleaning filenames to match popmap later on
for f in $out_dir/*.tsv*;
do
    mv $f $(echo "$f" | sed 's/-BWA-uniq.sorted//');
done;
