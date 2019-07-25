#!/bin/bash
# Downloads fastq files matching a range of identifiers on SRA
# Requires the fasterq-dump from NCBI
# cmdoret, 20190509

SRA_START=$1
SRA_END=$2
OUT_DIR=$3

# Extract longest common prefix between identifiers
SRA_PREFIX=$(grep -zPo '(.*).*\n\K\1' <<< "$SRA_START"$'\n'"$SRA_END")
SRA_NUM_START=${SRA_START#$SRA_PREFIX}
SRA_NUM_END=${SRA_END#$SRA_PREFIX}

# Iterate over range of variable ids and prepend common prefix
for ((sra_id=SRA_NUM_START;sra_id<=SRA_NUM_END;sra_id++)); do
		# Download each fastq
		fasterq-dump -e 4 "${SRA_PREFIX}"$(printf "%02" "$sra_id") -o $OUT_DIR
done
