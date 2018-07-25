# Using busco to assess completeness of the genome
# Cyril Matthey-Doret
# 17.02.2018

# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -r -d [-l] [-h]
   -r   reference genome to be assessed
   -d   BUSCO directory
   -l   run locally instead of using LSF platform.
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":r:ld:h" opt; do
   case $opt in
   r )  REF=${OPTARG} ;;
   d )  BUSCO_DIR=${OPTARG} ;;
   l )  local=yes;;
   h )  usage ;;
   \?)  usage ;;
   esac
done

if [ -z ${local+x} ]; then run_fun="bsub"; else run_fun="bash";fi

eval $run_fun << JOB

#!bin/bash
#BSUB -q normal
#BSUB -J BUSCO
#BSUB -e data/logs/busco.err
#BSUB -o data/logs/busco.out
#BSUB -M 32000000
#BSUB -R "rusage[mem=32000]"
#BSUB -n 30
#BSUB -R "span[ptile=30]"

source src/misc/dependencies.sh

# Prepare output folder
rm -rf $BUSCO_DIR
mkdir -p $BUSCO_DIR

# Retrieve database
wget -O $BUSCO_DIR/hymenoptera_odb9.tar.gz http://busco.ezlab.org/datasets/hymenoptera_odb9.tar.gz
tar -xzf $BUSCO_DIR/hymenoptera_odb9.tar.gz -C $BUSCO_DIR/

# run busco
run_BUSCO.py -i $REF -o lf_busco \
             -l $BUSCO_DIR/hymenoptera_odb9 \
             -m geno -c 30 -sp nasonia -f
mv run_lf_busco $BUSCO_DIR/lf_busco
JOB

