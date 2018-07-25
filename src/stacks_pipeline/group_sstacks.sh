
# This script is used for copying all snps, tags and alleles files from valid pstacks samples to sstacks folder
# It takes 4 arguments: the paths of P, C and S stacks folders, respectively, and whether (T) or not (F) it should 
# separate STACKS files into subfolders by family.
# Cyril Matthey-Doret
# 08.08.2017

PSTACK=$1
CSTACK=$2
SSTACK=$3

if [[ $# -eq 0 ]] ; then
    echo 'I need 3 arguments ! The path to pstacks, cstacks and sstacks output files, respectively. Exiting.'
    exit 0
fi

rm -f $SSTACK/batch*  # Remove files from previous runs
for f in $SSTACK/*;   # Iterating over sstacks files
do
    tf=$(basename ${f%%.matches*});   # Extracting individual names
    cp $PSTACK/$tf* $SSTACK   # Copying pstacks files corresponding to names to sstacks folder
done;
cp $CSTACK/* $SSTACK;  # Copying catalogue as well

