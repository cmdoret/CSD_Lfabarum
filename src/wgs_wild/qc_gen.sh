# Generate QC reports for both fastq files and Mapping.
# Needs to be run after the mapping step.
# Intended for the analysis of wild samples WGS data.
# Cyril Matthey-Doret
# 23.02.2018


function usage () {
  cat <<EOF
Usage: `basename $0` -d work_dir -r ref [-l] [-h]
  -d working directory. Must contain "raw" and "mapped" folders with input files.
  -r path to reference genome
  -l local run. If specified, will not use LSF bsub command
  -h displays this help
EOF
  exit 0
}

# Parse CL arguments
while getopts ":d:r:lh" opt; do
    case $opt in
        d ) wgs="${OPTARG}";;
        r ) index="${OPTARG}";;
        l ) local=yes;;
        h ) usage;;
        \?) usage;;
    esac
done

# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi

source src/misc/jobs_manager.sh
source src/misc/dependencies.sh
MAXPROC=15

# input and output folders
in_dir="${wgs}/raw/"
map_dir="${wgs}/mapped/"
out_dir="${wgs}/qc_output/"

# Resetting temporary working directory
tmp_dir="${wgs}/tmp/"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"

logs="${wgs}/log/"
mkdir -p "$logs"


# List of sample names
samples=( $(find "$in_dir" | grep "fastq" | sed 's/.*\(HYI-[0-9]*\)_R.*/\1/g' | sort | uniq) )

for sample in ${samples[@]};
do

if [ -z ${local+x} ];then bmonitor "WGSQC" $MAXPROC; fi

# Submit the QC of each sample as an independent job
eval $run_fun <<EOF
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o ${logs}/QC-${sample}-OUT.txt
#BSUB -e ${logs}/QC-${sample}-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSQC-${sample}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -q normal
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000

# Loading softwares
source src/misc/dependencies.sh

sample_fq="${tmp_dir}/${sample}.fastq.gz"

# Merge all fastq files for each sample (both lanes, reverse and forward)
find "${in_dir}" -name "*${sample}*" -type f | sort | xargs cat > "\$sample_fq"

# Fastqc on each sample
fastqc \$sample_fq -o $out_folder

# mapping stats with samtools stats
samtools stats -r $index "${map_dir}/${sample}.fixed.csort.bam" > $out_folder/${sample}_stat.txt

EOF
done

if [ -z ${local+x} ];then bmonitor "WGSQC" 0; fi

# MultiQC for all samples: Incorporating samtools picard and fastqc output
multiqc ${out_folder}/*  -o $out_folder

# remove temporary fq
rm -rf "$tmp_dir"
