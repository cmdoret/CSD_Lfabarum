#!/bin/bash
# Transforms basepairs into cM in the case_control_all.tsv file 
# from the association mapping using the output of bp2cm.sh
# Assumes both files are sorted by ascending Chr, BP
# TODO: distance between pairs of SNPs needs to be ratio*D
# Where D is the distance between a SNP and end of interval.
# If SNPs on different intervals, compute D per interval and sum them
# TODO: increment cumulative sum after every pair of SNP to get actual cM unit.
source src/misc/jobs_manager.sh
in_hits="$1"
bp2cm="$2"
out_hits="$3"

# Sort hits by position
hits="$in_hits.sorted"
head -n 1 "$in_hits" > "$hits"
sort -k2,2 -k3,3n <(tail -n+2 "$in_hits") >> "$hits"

# Cumulative cM value from previous markers
cum_cM=0
#Start at line 2 to skip header
n_interv=1
echo -e "$(head -n 1 "$hits")\tcM" > $out_hits

prev_chr=$(head -n 1 "$hits" | cut -f2)
# Loop over association mapping hits
while read -a hit; do
    ((i++))
    #prettyload "$i" $(wc -l "$hits" | awk '{print $1}')
    chr=${hit[1]}
    bp=${hit[2]}

    # Safety to reset the cM count on new chrom
    # works even if hit and interval changed chrom at the same time
    if [ $chr != $prev_chr ];then
        cum_cM=0
    fi
    # Loop over linkage intervals 
    # Stop reading if interval starts after hit on same chrom
    while read -a interv && ( ( [ ${interv[0]} == $chr ] && \
                                [ ${interv[1]} -le $bp ] ) || \
                              [ ${interv[0]} \< $chr ] ); do

        # If on different chrom, shift interval and reset cM value
        if [ ${interv[0]} \< $chr ];then
            cum_cM=0
            ((n_interv++))
            continue
        fi

        # If hit is inside segment (i.e. before its end)
        if [ $bp -le ${interv[2]} ]; then
            # Position of SNP from start of interval
            from_start=$((bp - interv[1]))
            curr_cM=$(echo "$from_start * ${interv[3]} + $cum_cM" | bc -l)
            echo -e "${hit[*]}\t$curr_cM" >> "$out_hits"
        else
            # SNP not found, next read will not use this interval
            ((n_interv++))
            # Computing total cM in interval (cM/BP) * (end-start) 
            int_cM=$(echo "${interv[3]} * (${interv[2]} - ${interv[1]})")
            cum_cM=$(echo "$cum_cM + $int_cM" | bc -l)
        fi
        echo "$cum_cM"

    # Only check intervals from n_inter to the end
    done < <(sed -n "${n_interv},\$p" $bp2cm)

    prev_chr=$chr
done < <( tail -n+2 $hits)

rm "$hits"
sed 's/ /	/g' "$out_hits" > "$out_hits.tmp" && mv "$out_hits.tmp" "$out_hits"
