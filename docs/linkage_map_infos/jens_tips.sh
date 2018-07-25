*****
clean run-through:

-set one-liner to remove sequence after MseI Adapter if run-through

cat file.fq | paste - - - - | tr '\t' 'z' | sed 's/TTAA.*z+z/TTAAz+z/' | tr 'z' '\t' | awk '{x = length($2); $4 = substr($4,0,x)}1' | tr ' ' '\n' > file.clip.fq







************************************
************************************
CORRECT WAY TO GET CLEAN 012 FILE:

-bowtie2 with RG flag
for f in *.fq.gz; do bowtie2 --end-to-end --very-sensitive --rg-id $f --threads 4 -x ../mapping/Lys -U $f -S $f.sam; done

-renaming RG info in read header and file header
(cat 001.fq.gz.sam | sed 's/001.fq.gz/001/g' | sed 's/@RG\tID:001/@RG\tID:001\tSM:001\tLB:1/g' > 001.sam)
(all in once)
ls *.sam > filelist
while read line; do IFS='.' read -a var <<< $line; cat ${var[0]}.${var[1]}.${var[2]}.${var[3]} | sed "s/${var[0]}.${var[1]}.${var[2]}/${var[0]}/g" | sed "s/@RG\tID:${var[0]}/@RG\tID:${var[0]}\tSM:${var[0]}\tLB:1/g" > ${var[0]}.sam; done < filelist

-convert to bam
for f in *.sam; do samtools view -bS -@ 24 $f -o $f.bam; done
-sort
for f in *.sam.bam; do samtools sort -@ 24 $f $f.sorted; done



-old pipe samtools (0.1.19) and bcftools:
time ~/Software/samtools-0.1.19/samtools mpileup -u -D -q 20 -f ../genome/Lys_final_assembly.fasta -b files.txt | ~/Software/samtools-0.1.19/bcftools/bcftools view -c -g -v -s names_ploidy.txt - > var.vcf
(on cluster:)
time samtools mpileup -u -D -q 20 -f ../canu2_low_corrRE.contigs.fasta -b files.txt | bcftools view -c -g -v -s names_ploidy - > var.vcf

files= list of input files
(ls *.sam.bam.sorted.bam > files.txt)
names_ploidy.txt= file with sample\t1 (haploid)
(ls *.sam | cut -d "." -f 1 | awk '{print$1,"\t1"}' | sed 's/ //g' > names_ploidy)


VCFTOOLS:
~/Software/vcftools_0.1.13/bin/vcftools --vcf test_new.vcf --minGQ 20 --max-missing 0.8 --maf 0.15 --thin --recode --recode-INFO-all --012 --out test_new_filtered


=> WORKS NOW!!! (RG FLAG!)




-> done with 90 biggers (>30MB) samples(individuals)
-> After filtering, kept 1319 out of a possible 868513 Sites



******************
preparing MSTmap:

- put "i" in front of sample individuals
cat Lys_var_filtered.012 | sed 's/^/i/'

-transpose rows&columns
~/Dropbox/scripts/transpose.sh test | sed 's/ /\t/g' > Lys_var_filtered_012.tp

-add loci
cat Lys_var_filtered.012.pos | sed 's/\t/_/g' | sed '1s/^/locus_name\n/' > loci
paste loci Lys_var_filtered_012.tp > Lys_MSTmap_in

-change -1,0,2 to U,A,B
cat Lys_MSTmap_in | sed 's/\t-1/\tU/g' | sed 's/\t0/\tA/g' | sed 's/\t2/\tB/g' > Lys_MSTmap.A

-double and reverse phase
cat Lys_MSTmap_in | sed 's/\t-1/\tU/g' | sed 's/\t0/\tB/g' | sed 's/\t2/\tA/g' > Lys_MSTmap.B

-Reversed Phase file and delete first ID line (locus ID will be A/B switched in name)
cat Lys_MSTmap.B | sed '1d' > Lys_MSTmap.B2

(test: apply tag "RP_" to locus ID:)
(cat Lys_MSTmap.B | sed '1d' | sed 's/^/RP_/' > Lys_MSTmap.B3)

-join both
cat Lys_MSTmap.A Lys_MSTmap.B2 > tmp

-add header (after setting parameters in txt file)
cat input tmp > Lys_MSTmap.inA



-input parameters in extra file (input):

population_type DH
population_name Lys
distance_function kosambi
cut_off_p_value 1e-6
no_map_dist 15.0
no_map_size 2 #set to 0 to disable this
missing_threshold 1.00
estimation_before_clustering no
detect_bad_data yes
objective_function COUNT
number_of_loci 2209 #(double number of loci because doubled, .i.e: wc -l tmp (-1 header line))
number_of_individual 124 #(number of samples)




*********
running MSTMap

~/Software/MSTMap/MSTMap.exe Lys_MSTmap.inA Lys_MSTmap_p1e-6.out

cut_off_p_value 1e-6


-> gives 5 linkage groups (= chr number Nasonia. L. fabarus might have 6 (paper))




Lys_MSTmap_p1e-6.out

- removing double Group (A/B in ID name switched)

-> Lys_curated_p1e-6.LG





search for microsattelites:

grep -c -i 'gacgctcagt' *.* | awk -F ':' '{if ($2==1) print$0}'

(mabye reverse complement)




************************************
************************************





parameters

samtools0.19
-u uncompressed bcf (Compute genotype likelihoods)
-f faidx index ref file
-D Output per-sample read depth 
-q Minimum mapping quality for an alignment to be used

bcftools0.19
-c Call variants using Bayesian inference
-v output variant sites only
-g Call per-sample genotypes at variant sites (force -c)
-s List of samples to use. The first column in the input gives the sample names and the second gives the ploidy, which can only be 1 or 2 (tab delimited). In the output, the ordering of samples will be identical to the one in FILE.



samtools 1.12
-u uncompressed
-g Compute genotype likelihood bcf
-f faidx index ref file
-t DP number of high quality bases (per sample read depth) 
-q Minimum mapping quality for an alignment to be used
-b list of input bam files

bcftools 1.12
-v output variant sites only 
-m multiallelic and rare-variant calling
-O v output type uncompressed vcf
-S file of sample names (second (tab del) column indicating ploidy (0, 1 or 2))
-o output



VCF tools parameters

--minGQ <float>
Exclude all genotypes with a quality below the threshold specified. This option requires that the "GQ" FORMAT tag is specified for all sites.

--max-missing <float>
Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).

--maf <float>
--max-maf <float>
Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value and less than or equal to the "--max-maf" value. One of these options may be used without the other. Allele frequency is defined as the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that site.

--thin <integer>
Thin sites so that no two sites are within the specified distance from one another.

--012
This option outputs the genotypes as a large matrix. Three files are produced. The first, with suffix ".012", contains the genotypes of each individual on a separate line. Genotypes are represented as 0, 1 and 2, where the number represent that number of non-reference alleles. Missing genotypes are represented by -1. The second file, with suffix ".012.indv" details the individuals included in the main file. The third file, with suffix ".012.pos" details the site locations included in the main file.

--recode
--recode-bcf
These options are used to generate a new file in either VCF or BCF from the input VCF or BCF file after applying the filtering options specified by the user. The output file has the suffix ".recode.vcf" or ".recode.bcf". By default, the INFO fields are removed from the output file, as the INFO values may be invalidated by the recoding (e.g. the total depth may need to be recalculated if individuals are removed). This behavior may be overriden by the following options. By default, BCF files are written out as BGZF compressed files.

--recode-INFO-all
These options can be used with the above recode options to define an INFO key name to keep in the output file. This option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file.




