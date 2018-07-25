#!/bin/bash
# This script handles the mapping of whole genome
# sequencing illumina reads from wild samples.
# Cyril Matthey-Doret
# 11.10.2017


function usage () {
  cat <<EOF
Usage: `basename $0` -d work_dir -r ref [-l] [-h]
  -d working directory. Must contain a "raw" directory with input reads.
  -r path to reference genome
  -l local run. If specified, will not use LSF bsub command
  -h displays this help
EOF
  exit 0
}

while getopts ":d:r:lh" opt; do
    case $opt in
        d ) wgs="${OPTARG}";;
        r ) ref="${OPTARG}";;
        l ) local=yes;;
        h ) usage;;
        \?) usage;;
    esac
done

# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi


# #threads used when mapping
threads=8
# Max number of parallel jobs
MAXPROC=30

# Allows to cap # parallel jobs and wait for their execution
source src/misc/jobs_manager.sh

in_dir="${wgs}/raw/"

# Each sample is split into forward and reverse on 2 different lanes
# Building array of sample names from raw reads files
samples=( $(find "$in_dir" \
            | grep "fastq" \
            | sed 's/.*\(HYI-[0-9]*\)_R.*/\1/g' \
            | sort \
            | uniq) )

# Reinitializing folders
merged_dir="${wgs}/merged/"
tmp_dir="${wgs}/tmp/"
map_dir="${wgs}/mapped/"
logs="${wgs}/log/"
for dir in "$merged_dir" "$tmp_dir" "$map_dir" "$logs"
do
  rm -rf "$dir"
  mkdir -p "$dir"
done


for sample in ${samples[@]};
do

# Hang script if too many parallel jobs running
if [ -z ${local+x} ];then bmonitor "WGSBWA" $MAXPROC;fi

# Submit the mapping of each sample as an independent job
eval $run_fun << EOF
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o ${logs}/${sample}-OUT.txt
#BSUB -e ${logs}/${sample}-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSBWA-${sample}
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -q normal
#BSUB -R "rusage[mem=64000]"
#BSUB -M 64000000

# Loading softwares
source src/misc/dependencies.sh

forward="${wgs}/merged/${sample}_R1.fastq.gz"
reverse="${wgs}/merged/${sample}_R2.fastq.gz"

# Merge lanes, keep forward and reverse separated
find "${in_dir}" -name "*${sample}*R1*" -type f | \
  sort | xargs cat > "\$forward"

find "${in_dir}" -name "*${sample}*R2*" -type f | \
  sort | xargs cat > "\$reverse"

# Trimmed reads filenames

# Paired:
trimF=\$(echo \$forward | sed 's%/\([^/]*\)$%/trim.\1%')
trimR=\$(echo \$reverse | sed 's%/\([^/]*\)$%/trim.\1%')

# Unpaired:
trimFU=\$(echo \$trimF | sed 's/R1.fastq.gz$/R1U.fastq.gz/')
trimRU=\$(echo \$trimR | sed 's/R2.fastq.gz$/R2U.fastq.gz/')

# Trim low quality ends
trimmomatic PE \$forward \$reverse \$trimF \$trimFU \$trimR \$trimRU \
            -trimlog ${logs}/${sample}-trim.log \
            LEADING:20 TRAILING:20

# Mapping paired ends reads
bwa mem -t $threads $ref \$trimF \$trimR > "${map_dir}/${sample}.sam"

# Convert SAM files to BAM
samtools view -@ $threads -b -o "${map_dir}/${sample}.bam" "${map_dir}/${sample}.sam"

# Sort alignments by read name
samtools sort -@ $threads -n "${map_dir}/${sample}.bam" -o "${map_dir}/${sample}.nsort.bam"

# Fix mate information (adds ms and MC tags for markdup)
samtools fixmate "${map_dir}/${sample}.nsort.bam" "${map_dir}/${sample}.fixed.bam"

# Sort alignments by leftmost coordinate
samtools sort -@ $threads "${map_dir}/${sample}.fixed.bam" \
              -o "${map_dir}/${sample}.fixed.csort.bam"
# Index BAM files
samtools index "${map_dir}/${sample}.fixed.csort.bam"

# Remove sam and unsorted/temporary bam files
rm -v "${map_dir}/${sample}.sam" \
      "${map_dir}/${sample}.bam" \
      "${map_dir}/${sample}*nsort*"

EOF
done

# Hang script while there are still mapping jobs running on cluster
if [ -z ${local+x} ];then bmonitor "WGSBWA" 0; fi
