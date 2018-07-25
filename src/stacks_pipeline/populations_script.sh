#!/bin/bash

# This script uses the populations component from the STACKS suite to compute various populations statistics the genome of individuals.
# It needs a popmap file as specified in the STACKS documentation: catchenlab.life.illinois.edu/stacks/comp/populations.php
# It requires stacks files for each individual in the same folder.
# Cyril Matthey-Doret
# 24.04.2017

source src/misc/jobs_manager.sh

while [[ "$#" > 1 ]]; do case $1 in
    # Path to popmap
    --pom) pom="$2";;
    # Path to Sstacks files
    --sst) sst="$2";;
    # Path to "threshold" file with ploidy information
    --thresh) thresh="$2";;
    # Minimum locus depth required
    --md) D="$2";;
    # Minimum proportion of individuals with locus
    --r) R="$2";;
    # Output populations folder
    --out) out_dir="$2";;
    # Location of log files
    --log) logs="$2";;
    # Blacklist of haploid loci
    --blacklist) black="$2";;
    # optional: run locally
    --local) local=yes;;
    *) break;;
  esac; shift; shift
done

# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi

# If popmap not provided, use default path
if [ -z ${local+x} ];then pom='data/popmap.tsv';fi

# Declare bl variable only if blacklist file exists at input path
if [ -f $black ]; then bl=$black; fi

    if [ -f $thresh ]  # If ploidy information already available
        then
            haplo=$(awk '$11 ~ "H" {print $1}' $thresh)  # list haploid males
            for indv in $haplo;
            do
                rm $sst/$indv*
                # Remove all files related to haploids to exclude from populations run
            done
            echo "Populations running on all individuals, excluding $(echo $haplo | wc -w) haploid males."
        fi
    mkdir -p $out_dir/  # Prepare output folder
    eval $run_fun <<POP
    #!/bin/bash
    #BSUB -J POP
    #BSUB -q normal
    #BSUB -e $logs/populations/POP_STDERR.err
    #BSUB -o $logs/populations/POP_STDOUT.out
    #BSUB -M 16000000
    #BSUB -R "rusage[mem=16000]"
    #BSUB -n 3
    #BSUB -R "span[ptile=3]"

    # Loading softwares
    source src/misc/dependencies.sh
    populations -P $sst -M $pom -p 2 -m $D -b 1 -r $R -k -f p_value \
    -t 3 --verbose --fstats --vcf --genomic --min_maf 0.1 --renz ecoRI ${bl:+-B "$bl"}

    mv $sst/batch* $out_dir/
    # Moving all populations output file from sstacks folder to populations folder
POP

# Only exit when populations jobs are finished
if [ -z ${local+x} ];then bmonitor POP 0; fi
