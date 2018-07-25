#!/bin/bash
# Genome-wide variant calling and formatting of output for downstream analyses.
# Cyril Matthey-Doret
# 02.01.2018

### HANDLE CL ARGUMENTS ##

# Help message
function usage () {
  cat <<EOF
Usage: `basename $0` -d work_dir -r ref [-l] [-h]
  -d working directory. Must contain a "mapped" directory with input bam files.
  -r path to reference genome
  -l local run. If specified, will not use LSF bsub command
  -h displays this help
EOF
  exit 0
}

# Parsing CL arguments
while getopts ":d:r:lh" opt; do
    case $opt in
        d ) WGS="${OPTARG}";;
        r ) REF="${OPTARG}";;
        l ) local=yes;;
        h ) usage;;
        \?) usage;;
    esac
done

# If on cluster, use bsub to submit jobs, otherwise run directly
if [ -z ${local+x} ];then run_fun="bsub";else run_fun="bash";fi


### SETUP ###
# Load environment and generate temporary files

# Clean and generate output folder
snps="$WGS/variant/"
mkdir -p "$snps"
rm -rf "${snps}/*"

logs="$WGS/logs/"
mkdir -p "$logs"


# Loading softwares
source src/misc/dependencies.sh
source src/misc/jobs_manager.sh

# Storing name and length of each contig in a temporary file
cat $REF | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; }
                $0 !~ ">" {c+=length($0);}
                END { print c; }' > tig_sizes.tmp
# Removing empty lines
sed '/^\s*$/d' tig_sizes.tmp > ts.bak && mv ts.bak tig_sizes.tmp

# Progress tracker (for loading bar)
tot_nucl=$(awk '{ sum += $2 } END { print sum }' tig_sizes.tmp)
done_nucl=0

# Chunk size (BP) to parallelize SNP calling by region on cluster
if [ -z ${local+x} ]
then
  CHUNK_SIZE=100000
else
  # if not on cluster, don't chunk (i.e. 1 large chunk)
  CHUNK_SIZE=$((10*$tot_nucl))
fi

# List of BAM files to read
find  "${WGS}/mapped/" -name "*fixed.csort.bam" -type f | \
	sed 's-//-/-g' > "${WGS}/bam_list.txt"

# List of samples filenames with corresponding ploidies
awk -v m="$WGS/mapped/" \
  'BEGIN{FS="\t"}
  {
    if($2 == "F") {p = 2;}
    else {p = 1}
    {print m$1".fixed.csort.bam",p}
  }' ${WGS}/wgs_samples.tsv | sed 's-//-/-g' > "${WGS}/spl_list.txt"

### RUN SNP CALLING ###
# SNP calling is performed on regions  of CHUNKSIZE bp
# in parallel to speed up the operation

# Iterating over contigs
while read -a line
do

  tig_ID=${line[0]}
  tig_len=${line[1]}
  start=1;end=0

  # Splitting current contig into regions
  while [ $end -lt $tig_len ]
  do

    # end takes smallest value between start+100kb and end of contig
    end=$(($tig_len<$end+$CHUNK_SIZE?$tig_len:$end+$CHUNK_SIZE))
    region=$tig_ID":"$start"-"$end

    # Track progress with a loading bar
    let "done_nucl += (end - start)"
    prettyload $done_nucl $tot_nucl

# Do not submit more than 200 jobs at once on cluster
if [ -z ${local+x} ];then
  bmonitor WGSSNP 200;fi

# 1 job per region
eval $run_fun << VAR > /dev/null
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o $logs/wgs_snp-OUT.txt
#BSUB -e $logs/wgs_snp-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSSNP-$region
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -q priority
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000

# Loading softwares
source src/misc/dependencies.sh

# SNP calling in given region for all samples (keeps only variant sites)
# File names are provided with ploidy by piping the sample list through awk
# Skipping indels because we are only interested in SNPs currently
bcftools mpileup -Ou -I \
                 -f "$REF" \
                 -b "${WGS}/bam_list.txt" \
                 -r "$region" | \
  bcftools call -mO z -M \
                --skip-variants indels \
                --samples-file "${WGS}/spl_list.txt" \
                -o "${snps}/$region.tmp.vcf.gz"

# Index vcf file
tabix -p vcf "${snps}/$region.tmp.vcf.gz"

VAR

    let "start += $CHUNK_SIZE"
done
done < tig_sizes.tmp

rm tig_sizes.tmp

# Wait for all SNP calling jobs to finish before resuming
if [ -z ${local+x} ];then
  bmonitor "WGSSNP" 0;fi

### PROCESS SNPS ###

# Clean and generate output folder
stats=${WGS}/stats
rm -rf $stats
mkdir -p $stats

eval $run_fun <<PROCSNP
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o $logs/proc_snp-OUT.txt
#BSUB -e $logs/proc_snp-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J PROCSNP
#BSUB -n 16
#BSUB -R "span[ptile=1]"
#BSUB -q normal
#BSUB -R "rusage[mem=14000]"
#BSUB -M 14000000

# Loading softwares
source src/misc/dependencies.sh

# format SNPs for downstream analyses
## 1: Concatenating all regions VCF files into a large one
find ${snps}/ -name '*tmp.vcf.gz' -print0 | \
  xargs -0 vcf-concat > ${snps}/wild.vcf

## 2: Parallel sorting of concatenated VCF file
vcf-sort -p 16 < ${snps}/wild.vcf > ${snps}/wild.sorted.vcf

## 3: Generate SNP matrix from VCF
bcftools query ${snps}/wild.sorted.vcf -f '%CHROM\t%POS[\t%GT]\n' | \
        grep -v "^##" > ${snps}/wild.matrix.txt

## 4: separate haplotypes of diploid samples into 2 columns
## Keeping only sites anchored to a chromosome
cut -f1-2 ${snps}/wild.matrix.txt > left_col.txt
cut -f3- ${snps}/wild.matrix.txt | sed 's#/#\t#g' > right_col.txt
paste left_col.txt right_col.txt | sed -n '/^chr/p' > ${snps}/hap.wild.matrix.txt
rm left_col.txt right_col.txt

## 5: Fill missing positions (if any) to make downstream analyses more convenient
python2 $(dirname $0)/fill_pos_matrix.py ${snps}/hap.wild.matrix.txt > mat.tmp && \
  mv mat.tmp ${snps}/hap.wild.matrix.txt

## 6: Generate second matrix with nucleotides, repeating steps 3-5 with a different
bcftools query ${snps}/wild.sorted.vcf -f '%CHROM\t%POS[\t%TGT]\n' | \ 
        grep -v "^##" > ${snps}/nuc.matrix.txt

cut -f1-2 ${snps}/nuc.matrix.txt > left_col.txt
cut -f3- ${snps}/nuc.matrix.txt | sed 's#/#\t#g' > right_col.txt
paste left_col.txt right_col.txt | sed -n '/^chr/p' > ${snps}/hap.nuc.matrix.txt

rm left_col.txt right_col.txt
python2 $(dirname $0)/fill_pos_matrix.py ${snps}/hap.nuc.matrix.txt > mat.tmp && \
  mv mat.tmp ${snps}/hap.nuc.matrix.txt
PROCSNP

if [ -z ${local+x} ];then
  bmonitor PROCSNP 0;fi

# Removing all regions VCF
#rm ${snps}/*.tmp.vcf.gz* ${WGS}/bam_list.txt ${WGS}/spl_list.txt
