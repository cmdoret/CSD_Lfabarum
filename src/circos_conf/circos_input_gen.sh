# This script generates the input files required to run circos.
# Cyril Matthey-Doret
# 21.11.2017

#---SET UP DIRECTORY---#
# Reference genome
ref="data/ref_genome/ordered_genome/merged.fasta"
# Correspondance file between old and new (anchored) contigs
tig="data/annotations/corresp_gff.csv"
# All SNPs from the association mapping
gwas="data/assoc_mapping/case_control/case_control_all.tsv"
# Significant SNPs from the association mapping
gwas_hit="data/assoc_mapping/case_control/case_control_hits.tsv"
# Collinearity blocks from MCScanX
mcsx="data/homology/MCScanX/input/MCScanX_in.collinearity"
# GFF with transcripts coordinates
gff="data/homology/MCScanX/input/MCScanX_genes_conv.gff"
cir_dir="data/circos/"
rm -rf "$cir_dir"
mkdir -p "$cir_dir"
# power of 10 to consider hits significant
sigpow=5

#---FORMAT DATA FILES---#

# Karyotype file:
# Parsing FASTA reference to karyotype format: chr - name label start end color
awk '$0 ~ ">" {
         print c; c=0;printf substr($0,2,100) "\t"; }
     $0 !~ ">" {
         c+=length($0);}
       END { print c; }' $ref |
  grep "chr" |
  sed 's/chr/lf/' |
  awk 'BEGIN{
         OFS=" ";N=0}
       {
         N+=1;print "chr","-",$1,N,0,$2,"chr"N}' > "$cir_dir/karyotype.lf.txt"

# Contig boundaries
tail -n +2 $tig |
  awk 'BEGIN{
         FS=",";OFS=" "}
       $2 ~ "chr" {
         if ( $5 ~ "rev" ){print $2,$3-$4,$3}
       else {print $2,$3,$3+$4}
       }' |
  sed 's/chr/lf/' |
  sort -k1,1 -k2,2n > "$cir_dir/tig_lim.lf.txt"

# Contigs with CSD hits
while read line
do
  # is the p-value significant ?
  pval=$(echo $line | awk '{print $NF}' | sed 's/[eE]+\{0,1\}/*10^/g')
  if (( $(echo "$pval  > $sigpow" | bc -l) ))
  then
    awk -v chrom="$(cut -d$' ' -f2 <(echo $line))" \
        -v pos="$(cut -d$' ' -f3 <(echo $line))" \
        '$1 == chrom {
          if ( $2 <= pos && $3 >= pos ) {print $0}}' "$cir_dir/tig_lim.lf.txt"
  fi
done < <( tail -n +2 "$gwas_hit" | sed 's/chr/lf/') |
  uniq > "$cir_dir/csd_tig.lf.txt"

# GWAS Manhattan plot
# parsing p-values of SNPS into: chr start end value
# where start=end (points)
tail -n +2 $gwas |
  awk 'BEGIN{IFS="\t";OFS=" "}
       {print $2,$3,$3,$13}' |
  grep "chr" |
  sed 's/chr/lf/' > "$cir_dir/gwas.lf.txt"

# Collinearity blocks
# Parsing MCScanX collinearity output into circos links
while read line
do
  awk -v g1="${line%% *}" -v g2="${line##* }" \
    'BEGIN{FS=" ";N1=0;N2=0}
     $2 ~ g1 {
       N1+=1;c1=$1;s1=$3;e1=$4}
     $2 ~ g2 {
       N2+=1;c2=$1;s2=$3;e2=$4}
     END{
       if (N1 != 1 || N2 != 1) {exit 1}
       else {print c1,s1,e1,c2,s2,e2}}' $gff
  # Convert gene IDs to coordinates
done < <(grep -v "^#" $mcsx | tr -d ' ' | tr '\t' ' ' | cut -d$' ' -f2,3) |
  sed 's/chr/lf/g' > "$cir_dir/mcsx.lf.txt"

# Alternative links: BLAST hits of transcripts
python2 src/circos_conf/blast_gene_circos.py
sed 's/chr/lf/g' "data/circos/blast.lf.txt" > "tmp" && \
  mv "tmp" "data/circos/blast.lf.txt"

python2 src/circos_conf/CSD_blast.py
# Computing homozygosity along chromosomes (centromere approximation)
Rscript src/circos_conf/centro_heatmap.R
sed 's/chr/lf/' "$cir_dir/centro.lf.txt" > "tmp_centro" && \
  mv "tmp_centro" "$cir_dir/centro.lf.txt"


#---GENERATE CONFIGURATION---#
# config file:
cat << CFG > $cir_dir/lf.main.conf
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<<include $cir_dir/lf.ideogram.conf>>
<<include $cir_dir/lf.ticks.conf>>
<<include $cir_dir/lf.mcscanx.conf>>
#<<include $cir_dir/lf.blast.conf>>
karyotype = "$cir_dir/karyotype.lf.txt"
chromosomes_units = 1000000
# Standard stuff:
<image>
<<include etc/image.conf>>
</image>
CFG

cat << IDEOGRAM > $cir_dir/lf.ideogram.conf
<ideogram>
<spacing>
# Spacing between ideograms. Suffix "r" denotes a relative value. It
# is relative to circle circumference (e.g. space is 0.5% of
# circumference).
default = 0.005r
#<pairwise lf1;lf2>
#spacing = 20r
#</pairwise>

</spacing>
radius           = 0.90r
thickness        = 60p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p

# LABELS
show_label       = yes
label_font       = default
label_radius     = 1r + 75p
label_size       = 50
label_parallel   = yes
label_format   = eval(sprintf("chr%s",var(label)))

</ideogram>

<highlights>

z          = 5

 <highlight>
 # Contig boundaries
 ideogram   = yes
 file       = "$cir_dir/tig_lim.lf.txt"
 fill_color = blue
 stroke_color = black
stroke_thickness = 2
 </highlight>

 <highlight>
 # CSD contig
 ideogram   = yes
 file       = "$cir_dir/csd_tig.lf.txt"
 z          = 5
 fill_color = red
 </highlight>
</highlights>

<plots>
<plot>
# CSD highlight
type = highlight
fill_color = red_a5
stroke_color = red_a5
file       = "$cir_dir/csd_tig.lf.txt"
r0   = 0.65r
r1   = 1r
z    = 10
</plot>
<<include $cir_dir/lf.centro.conf>>
<<include $cir_dir/lf.gwas.conf>>
</plots>

IDEOGRAM

cat << CENTRO > $cir_dir/lf.centro.conf
<plot>
show = yes
type = heatmap
file = $cir_dir/centro.lf.txt
r1= 0.70r
r0 = 0.65r
min = 0.3
max = 0.5
#color = blues-9-seq-rev
color = ylgnbu-6-seq

</plot>
CENTRO

cat << GWAS > $cir_dir/lf.gwas.conf
#GWAS p-values (scatter plot)
<plot>

show  = yes
type  = scatter

file  = $cir_dir/gwas.lf.txt
r1    = 1r - 55p
r0    = 0.70r + 10p
max   = 6.2
min   = 0.0
orientation = in

glyph            = circle
glyph_size       = 12
color            = vdgrey
stroke_color     = vdgrey
stroke_thickness = 1

<axes>
<axis>
color = vlgrey
thickness = 3
position = $sigpow
</axis>
</axes>

<backgrounds>
<background>
color = prgn-3-div-2
y0 = 0
</background>
</backgrounds>

<rules>
<rule>
condition    = var(value) > $sigpow
color        = dred
fill_color   = dred
glyph_size       = 18
</rule>
</rules>

</plot>

GWAS


cat << MCSCANX > $cir_dir/lf.mcscanx.conf
# Collinearity blocks
<links>

<link>
file          = "$cir_dir/mcsx.lf.txt"
radius        = 0.7r
bezier_radius = 0r
color         = lgrey_a4
thickness     = 3
</link>
</links>

MCSCANX

cat << BLAST > $cir_dir/lf.blast.conf
# BLAST hits outside CSD regions
<links>

<link>
file = $cir_dir/out_blast.lf.txt
radius        = 0.65r
bezier_radius = 0r
color         = lgrey_a4
thickness     = 2
<rules>
<rule>
condition     = var(intrachr)
show = no
</rule>
</rules>

</link>

<link>
# BLAST hits in CSD regions
file = $cir_dir/csd_blast.lf.txt
radius        = 0.65r
bezier_radius = 0r
color         = red_a3
thickness     = 8
z = 30

<rules>
<rule>
condition     = var(intrachr)
show = no
</rule>
</rules>

</link>
</links>

BLAST

# OPTIONAL
# coverage / FPKM ?
# allelic diversity
# transcripts



cat << TICKS > $cir_dir/lf.ticks.conf
# TICKS

show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = 1r
color            = black
thickness        = 2p

multiplier       = 1e-6

format           = %d

<tick>
show_label = yes
label_size = 25
spacing        = 5u
size           = 20p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>
</ticks>
TICKS

#---RUN CIRCOS---#
export PATH="~/Public/circos-0.69-6/bin/":$PATH
circos -conf "$cir_dir/lf.main.conf"
