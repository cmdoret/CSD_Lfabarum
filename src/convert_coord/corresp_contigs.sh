# This script allows to get a correspondance file of the coordinates of all
# contigs between 2 references. Both references must contain the same sequence,
# but contigs can be reorded, reverse-complemented and merged in the new
# reference.
# Cyril Matthey-Doret
# 31.10.2017

# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -O old_ref -N new_ref -c corresp [-l] [-h]
   -O   old reference
   -N   new reference
   -g   path to log file
   -c   output correspondance file
   -l   local run. If specified, will not use LSF bsub command
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":O:N:g:lc:h" opt; do
   case $opt in
   O )  OLD_REF=${OPTARG} ;;
   N )  NEW_REF=${OPTARG} ;;
   g )  LOGS=${OPTARG} ;;
   c )  CORRESP_GFF=${OPTARG};;
   l )  local=yes;;
   h )  usage ;;
   \?)  usage ;;
   esac
done

if [[ "x" == "x$OLD_REF" ]] || [[ "x" == "x$NEW_REF" ]] || [[ "x" == "x$CORRESP_GFF" ]];
then
  echo "Error: You must provide the path to both reference genomes as well as the output path for \
the correspondance file"
  usage
  exit 0
fi

tmp_dir=$(dirname $CORRESP_GFF)/tmp/
mkdir -p $tmp_dir


if [ -z ${local+x} ];then run_fun="bsub -I -tty";else run_fun="bash";fi

eval "$run_fun" <<CORR
#!/bin/bash
#BSUB -J corresp_GFF
#BSUB -q normal
#BSUB -e $LOGS/gff_corresp.err
#BSUB -o $LOGS/gff_corresp.out
#BSUB -M 32000000
#BSUB -R "rusage[mem=32000]"
#BSUB -n 36
#BSUB -R "span[ptile=36]"



# Get start position of all contigs in the ordered assembly and whether
# They have been reversed and/or complemented. Store infos in a file.

declare -i n=0
grep ">" $OLD_REF \
    | sed -n 's/>\([^ ]*\) .*/\1,0/p' \
    | python2 src/convert_coord/contig2chr.py \
                        --region_size 1000000000 \
                        --include_input \
                        $OLD_REF \
                        $NEW_REF \
    2> /dev/null \
    | cut -d ',' -f2 --complement \
    > $CORRESP_GFF

# insert header
sed -i "1s/^/tig,chr,start,region_size,transform\\n/" $CORRESP_GFF

CORR

rm -rf $tmp_dir
