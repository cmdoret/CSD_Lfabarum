# Parsing the output VCF file from populations runs in each family
# into summarized tables and convenient genotype matrices.
# Takes 2 argument: input path and whether populations was run separately
# for each family (F), or once for all individuals (T)
# Cyril Matthey-Doret
# 24.06.2017


if [[ $# -eq 0 ]]; then
    echo "No argument given, exiting now."
    echo "I need the folder where populations output files are stored and an output folder."
    exit 0
fi
declare -i increm=0
out_path=$2
rm -rf $out_path
mkdir -p $out_path

mkdir -p $out_path/all/

if [ -f $1/*.vcf ]
then
    vcftools --vcf $1/*.vcf --out $out_path/all/all --depth;
    # Mean read depth per individual
    vcftools --vcf $1/*.vcf --out $out_path/all/all --het;
    # Mean heterozygosity per individual
    vcftools --vcf $1/*.vcf --out $out_path/all/all --012;
    # Genotype matrix[SNP, individual]
    paste $out_path/all/all.het $out_path/all/all.idepth \
    | cut -f-5,8- > $out_path/summary_full.txt;
    # Putting together mean depth and mean heterozygosity
fi

# Removing empty lines from summary. bak file is used for BSD-sed compatibility
sed '/nan/d;' "$out_path/summary_full.txt" > "$out_path/summary_full.txt.tmp" && \
    mv "$out_path/summary_full.txt.tmp" "$out_path/summary_full.txt"
