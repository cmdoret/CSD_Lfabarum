#!/bin/bash
# This script will build a catalogue of loci using cstacks, given a set of pstacks output.
## It will exclude samples with very low number of reads from the catalogue (less than 10% of
## the mean number of reads across samples)
# Cyril Matthey-Doret
# 22.04.2017

source src/misc/jobs_manager.sh
# Default mismatches allowed between samples when building loci set to 1
LM=1
# n is the number of samples, tot will store total number of radtags
declare -i n=0 tot=0

while [[ "$#" > 1 ]]; do case $1 in
    # Path to Pstacks files
    --pst) pst="$2";;
    # Output path for catalog
    --cst) cst="$2";;
    # Loci mismatches allowed between alleles
    --lm) LM="$2";;
    # optional: run locally
    --local) local=yes;;
    *) break;;
  esac; shift; shift
done

# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi

for f in $pst/*tags*;  # Iterating over all samples
do
    tot+=$(gunzip -c $f | wc -l);  # adding number of radtags to total
    n+=1;  # incrementing number of samples by one at each iteration
done;


samp=""  # initiating list of "good" samples
for i in $pst/*tags*;  # Iterating over samples (again)
do
    if [ "$(gunzip -c $i | wc -l)" -gt $(($tot/($n*10))) ];
    then
        samp+="-s ${i%%.tags*} "  # only samples containing more than 10% of (arithmetic) mean radtags are used
    fi;
done;

eval $run_fun <<CST
#!/bin/bash

#BSUB -o %J_STDOUT.log
#BSUB -e %J_STDERR.log
#BSUB -J Cstacks_${j}
#BSUB -n 5
#BSUB -R "span[ptile=5]"
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10000000
#BSUB -q long

# Loading softwares
source src/misc/dependencies.sh

cstacks -b 1 -o $cst -n $LM -p 5 $samp;
CST

if [ -z ${local+x} ];then bmonitor Cstacks 0; fi
