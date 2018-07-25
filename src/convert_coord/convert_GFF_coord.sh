#!/bin/bash
# This script converts the coordinates in an original GFF file to coordinates
# corresponding to another assembly. It uses a correspondence file generated
# by corresp_contigs.sh to retrieve contig coordinates in the new assembly.
# Cyril Matthey-Doret
# 02.11.2017

# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -i input_gff -o output_gff (-c corresp_file | \
-O old -N new) [-l] [-h]
   -i   gff file to be converted
   -o   output converted gff file
   -c   csv file for correspondance between contigs
   -O   old reference, only needed if -c is not specified
   -N   new reference, only needed if -c is not specified
   -l   local run. If specified, will not use LSF bsub command
   -h   displays this help
EOF
   exit 0
}

# Parsing CL arguments
while getopts ":i:o:c:O:N:lh" opt; do
   case $opt in
   
   g )  LOGs=${OPTARG} ;;
   i )  GFF=${OPTARG} ;;
   o )  OUT_GFF=${OPTARG} ;;
   c )  CORRESP_GFF=${OPTARG};;
   O )  OLD_REF=${OPTARG};;
   N )  NEW_REF=${OPTARG};;
   l )  local=yes;;
   h )  usage ;;
   \?)  usage ;;
   esac
done

if [ "x" == "x$GFF" ] || [ "x" == "x$OUT_GFF" ];
then
  echo "Error: Input and output GFF files must be provided \
  along with either a contig correspondance file or both references."
  usage
  exit 0
fi

shift $(($OPTIND - 1))
# set command to be used depending if running locally or on LSF
if [ -z ${local+x} ];then run_fun="bsub -I -tty";else run_fun="bash";fi

# Generate correspondance file if not specified
if [ -z ${CORRESP_GFF+x} ];
then
  if [ -z "$OLD_REF" ] || [ -z "$NEW_REF" ];
  then
    echo "Error: If no correspondance file is provided, both the old and new \
    references are required."
    usage
    exit 0
  fi
  # Default path for corresp file
  CORRESP_GFF="data/annotations/corresp_gff.csv"
  bash corresp_contigs.sh -O "$OLD_REF" \
                          -N "$NEW_REF" \
                          -G "$LOGS"
                          -c "$CORRESP_GFF" ${local:+-l}
else if [ ! -f "$CORRESP_GFF" ];
  then
    echo "Correspondance file invalid. Exiting"
    exit 0
  fi
fi

# Run main code with bsub or directly, depending on -l flag
# Note: many variables are escaped to force remote expansion of heredoc
echo "Job submitted !"
eval $run_fun <<CONV_COORD

  #!/bin/bash

  #BSUB -J conv_GFF
  #BSUB -q normal
  #BSUB -e $LOGS/gff_conv.err
  #BSUB -o $LOGS/gff_conv.out
  #BSUB -M 16000000
  #BSUB -R "rusage[mem=16000]"
  #BSUB -n 28
  #BSUB -R "span[ptile=28]"

  source src/misc/jobs_manager.sh

  # Performance tuning parameters
  MAX_PROC=24;CHUNK_SIZE=200
  n_rec=\$(wc -l $GFF | awk '{print \$1}')
  n_chunk=\$((\$n_rec / \$CHUNK_SIZE + 1))

  # Clean temporary files
  tmpdir="\$(dirname $OUT_GFF)/tmp/"
  rm -rf \$tmpdir && mkdir -p \$tmpdir
  echo -n "" > $OUT_GFF

  for ((chunk=0;chunk < n_chunk; chunk++))
  do
    # Line boundaries of the chunk (1-indexed)
    c_start=\$(( \$chunk * \$CHUNK_SIZE + 1))
    c_end=\$((\$c_start + (\$CHUNK_SIZE - 1)))
    # Stop spawning subprocesses if too many running
    [ \$( jobs -p | wc -l ) -ge \$MAX_PROC ] && wait

    # spawning one parallel subprocess per chunk
    # each subprocess iterates over lines of its chunk
    (
    sed -n "\${c_start},\${c_end}p" $GFF | while read line
    do
      # Storing line of GFF file into bash array
      track=( \$line )
      # Getting line of corresp file for the contig on current line
      corresp=( \$(grep "^\${track[0]}" $CORRESP_GFF | sed 's/,/ /g') )
      # Replacing contig of current record in memory
      track[0]=\${corresp[1]}
      # Shift start and end if necessary and flip,
      # depending if contig was reversed or not
      if [[ \${corresp[4]} == *"rev"* ]]
      then
        # Reversed: corresp[2] will match the end of the contig
        # start -> contig_end-track_end, end -> contig_end-track_start
        start=\$((\${corresp[2]}-\${track[4]}))
        end=\$((\${corresp[2]}-\${track[3]}))
        track[3]=\$start
        track[4]=\$end
      else
        # not reversed: start shifted, end -> start + contig size
        let "track[3] += \${corresp[2]}"
        let "track[4] += \${corresp[2]}"
      fi
      # If complementary and not strand agnostic -> complement track
      if [[ \${corresp[4]} == *"comp"* && \${track[6]} != "." ]]
      then
        # Reverse strand
        if [[ \${track[6]} == "+" ]]
        then
          track[6]="-"
        else
          track[6]="+"
        fi
      fi    # redirect line to output gff (line >> file)
      # Write line to temporary file to avoid write conflicts
      echo "\${track[@]}" | tr ' ' \\\\t >> "\$tmpdir/chunk.\$chunk";
    done
    ) &
    # Displaying progress bar
    prettyload \$c_start \$n_rec
  done
  wait

  # Concatenate temporary files (this syntax allows to get over
  # cat's maximum number of argument)
  find \$tmpdir/ -name "chunk*" -type f -maxdepth 1 | \
      xargs cat > $OUT_GFF
  # Sort lines according to new coordinates
  sort -k1,1 -k4,4n -o $OUT_GFF $OUT_GFF
CONV_COORD
