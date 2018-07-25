#!/bin/bash
#BSUB -J radtags_a
#BSUB -o demulti-a-output.o
#BSUB -e demulti-a-error.e
#BSUB -n 4
#BSUB -u cmatthey@unil.ch
#BSUB -R "span[ptile=4]"
#BSUB -q normal
#BSUB -R "rusage[mem=4000]"
#BSUB -M 8000000

# This script processed the raw reads for an a-type library (adapter iA_06), using the process_radtags program from STACKS.
# Takes the library ID and optionally, the number of authorized mismatches in adaptor as an argument
# Cyril Matthey-Doret
# 15.03.2017

# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -l lib_path -a adapter_name -b barcodes.txt -o out_folder [-m mismatches] [-h]
Processes raw reads using the STACKS process_radtags module. Trims adapters and demultiplexes samples.
   -l   Path to the raw library to be processed
   -a   Name of the adapter sequence. iA_06 and iA_12 supported.
   -b   Barcodes file. List of tab separated barcodes and names, one sample per line.
   -o   Output folder where the processed reads will be written.
   -m   Number of mismatches allowed [2]
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":l:a:b:o:mh" opt;do
    case $opt in
        l ) LIB=${OPTARG} ;;
        a ) ADA=${OPTARG} ;;
        b ) BAR=${OPTARG} ;;
        o ) OUT=${OPTARG} ;;
        m ) MIS=${OPTARG} ;;
        h ) usage ;;
        \?) usage ;;
    esac
done

if [[ "x" == "x$LIB" ]] || [[ "x" == "x$BAR" ]] || [[ "x" == "x$OUT" ]]
then
        echo "You must specify the path to the raw library, barcodes and output folder."
        usage
fi

case "$ADA" in
        "iA_06" ) ADA="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG " ;;
        "iA_12" ) ADA="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG " ;;
        * ) echo "Adapter unsupported"; usage ;;
esac

if [[ "x" == "x$MIS" ]]
then
    MIS=2
fi


mkdir -p $OUT

source src/misc/dependencies.sh

process_radtags -r -q -c  \
                -e ecoRI  \
                --filter_illumina \
                -i gzfastq \
                -f "$LIB" \
                -b "$BAR" \
                -o "$OUT" \
                --adapter_1 "$ADA" \
                --adapter_mm $MIS

