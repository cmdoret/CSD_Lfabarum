
#!/bin/bash


# This script runs sstacks, matching samples against the catalogue.
# It will exclude samples with very low number of reads from the analysis (less than 10% of
# the mean number of reads across samples), as they were not used when building the catalogue
# Cyril Matthey-Doret
# 10.10.2017

source src/misc/jobs_manager.sh
declare -i n=0 tot=0 # Number of samples and total radtags

# parsing CL arguments
while [[ "$#" > 1 ]]; do case $1 in
    # Path to the folder containing input pstacks
    --in) pst="$2";;
    # Path to Cstacks catalog
    --cat) cst="$2";;
    # Path to output stacks matches files
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
rm -rf $logs/sstacks; mkdir -p $logs/sstacks

# Summing all Total reads and number of samples
for f in $pst/*tags*;
do
    tot+=$(gunzip -c $f | wc -l);
    n+=1;
done;

# Filtering out low quality samples
samp=""
for i in $pst/*tags*;
do
    if [ "$(gunzip -c $i | wc -l)" -gt $(($tot/($n*10))) ];
    # Only using samples containing more radtags than 10% of arithmetic mean over all samples
    then
        samp+="${i%%.tags*} "
    fi;
done;

# Writing scripts automatically
for i in $samp  # Writing a script automatically for each "good" sample
do
    # Do not queue more than 100 jobs at a time
    if [ -z ${local+x} ];then bmonitor SST 100; fi
    echo "Sample $(basename $i)"
    j=$(echo ${i##*/} | cut -f1 -d '.')

    eval $run_fun <<SST
    #!/bin/bash
    #BSUB -L /bin/bash
    #BSUB -o $logs/sstacks/SST_${j}_STDOUT.log
    #BSUB -e $logs/sstacks/SST_${j}_STDERR.log
    #BSUB -J SST${j}
    #BSUB -n 3
    #BSUB -M 2000000
    #BSUB -q priority

    # Loading softwares
    source src/misc/dependencies.sh
    sstacks -b 1 -c $cst/batch_1 -s $i -o $out_dir -p 3
SST
done

# Script waits for all sstacks jobs to finish
if [ -z ${local+x} ];then bmonitor SST 0; fi
