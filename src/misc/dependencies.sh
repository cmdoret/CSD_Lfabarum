# This script centralizes dependency path for the whole project
# Cyril Matthey-Doret
# 02.02.2018


# Replace paths appropriately
module add UHTS/Aligner/bwa/0.7.15 \
           UHTS/Analysis/samtools/1.4 \
           UHTS/Analysis/vcftools/0.1.15 \
           R/3.4.2 \
           Blast/ncbi-blast/2.6.0+ \
           UHTS/Analysis/BEDTools/2.26.0 \
           UHTS/Quality_control/fastqc/0.11.5 \
           UHTS/Analysis/MultiQC/1.3 \
           UHTS/Analysis/deepTools/2.5.4 \
           SequenceAnalysis/GenePrediction/augustus/3.2.3 \
           SequenceAnalysis/HMM-Profile/hmmer/3.1b2 \
           UHTS/Analysis/stacks/1.48 \
           UHTS/Assembler/cufflinks/2.2.1 \
           UHTS/Analysis/trimmomatic/0.36 \
           UHTS/Analysis/HTSlib/1.6 2> /dev/null || echo "Trying to load modules locally"

# Loaded from local install
SOFT="$HOME/softwares"
export PATH=$PATH:"$SOFT/MCScanX/"
export PATH=$PATH:"$SOFT/busco/3.0.2b/scripts/"
export PATH=$PATH:"$SOFT/cufflinks-2.2.1/"

# Necessary for VCFtools perl scripts
export PERL5LIB=/software/UHTS/Analysis/vcftools/0.1.15/lib64/perl5/site_perl/5.24.0/

# AUGUSTUS config required to use BUSCO
export AUGUSTUS_CONFIG_PATH="$SOFT/busco/3.0.2b/augustus_config/"

# Alias to java class for picard
function picard {
    java -jar -XX:MaxPermSize=32g "$SOFT/picard/build/libs/picard.jar" "$@"
}


LEPMAP3="$SOFT/lepmap3/bin"
