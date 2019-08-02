#!/bin/bash

GENOMESIZE=140734580
READLENGTH=93
REPORT="sample_coverage.tsv"

echo -e "sample\tcoverage" > "$REPORT"
n_samples=$(ls -1 data/mapped/aln-4/bam/*bam | wc -l)
declare -i n_done=0
for sample_bam in data/mapped/aln-4/bam/*bam; do
        echo -e "Processed $n_done / $n_samples"
        COVERAGE=$(samtools depth -a "$sample_bam" \
                | awk -v gsize="$GENOMESIZE" -v rlen="$READLENGTH" '{depth+=$3}END{print rlen*depth/gsize}')
        sample_name="$(basename $sample_bam)"
        echo -ne "${sample_name%%-BWA*}\t" >> "$REPORT"
        echo "$COVERAGE" >> "$REPORT"
        ((n_done++))
done
# Merge coverage file with sample metadata
join -j 1 \
     <( tail -n+2 sample_coverage.tsv | sort -k1,1) \
     <( tail -n +2 data/ploidy/thresholds/fixed.tsv | sort -k1,1)
