include config.mk

# LOCAL must be defined from command line if running analysis locally
# e.g: make collinearity LOCAL=yes
ifdef $(LOCAL)
		LOCAL='$(LOCAL)';
endif

##########################################
#### 1. RADseq processing with STACKS ####
##########################################

# reference-based STACKS pipeline for LSF cluster (default)
ref_lsf:
	make -f src/pipelines/Makeref.lsf

# reference based STACKS running locally (without LSF)
.PHONY : local
ref_local:
	make -f src/pipelines/Makeref.local

##############################
#### 2. PLOIDY SEPARATION ####
##############################

# This rule is used to split haploid and diploid males.
# This has already been done with the dataset provided
# under stringent parameters (D=25 in STACKS populations)
.PHONY : ploidy
ploidy:
	# Parsing VCF file from populations output
	mkdir -p $(DAT)/ploidy/thresholds
	bash $(MISC)/parse_VCF.sh $(POP) $(DAT)/ploidy/vcftools
	# Building list of haploid males
	python2 src/ploidy/haplo_males.py $(VCFSUM) \
			$(THRESH) \
			--sample_list $(DAT)/individuals.tsv \
			--ploidy_thresh $(HOM_PLOID)
	# Processing populations genomic output (excluding hom/missing SNPs in
	# mothers from their families)
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) \
			$(DAT)/ploidy/
	# Blacklisting loci that are heterozygous in haploid males
	python2 src/ploidy/blacklist_haploloci.py $(DAT)/ploidy/ \
			$(BLACK) \
			$(THRESH)
	mkdir -p data/SNP_lists

################################
#### 3. ASSOCIATION MAPPING ####
################################

# Association mapping
.PHONY : assoc_mapping
assoc_mapping : $(CENTRO) $(SIZES) $(ASSOC)/grouped_prophom.tsv
	# Creating folder to store association mapping results
	mkdir -p $(ASSOC)/plots
	mkdir -p $(ASSOC)/hits
	# Genome-wide association mapping
	mkdir -p $(ASSOC)/case_control
	Rscript $(ASSOC-SRC)/case_control.R \
			$(ASSOC)/grouped_outpool_keep_prophom.tsv \
			$(ASSOC)/case_control/ \
			$(SIZES)


# Centromere identification
$(CENTRO) :
	mkdir -p $(CENTRO)/plots
	# Processing "genomic" output from populations to get fixed and variant sites
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) \
			--pool_output
	# Inferring centromere position based on recombination raes along chromosomes
	Rscript $(ASSOC-SRC)/chrom_types.R $(ASSOC)/grouped_outpool_prophom.tsv \
			$(DAT)/individuals.tsv \
			$(CENTRO)

# Processing STACKS "genomic" output for association mapping
$(ASSOC)/grouped_prophom.tsv :
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) \
			--keep_all --pool_output


# Store sizes of all contigs in text file
$(SIZES):
	awk '$$0 ~ ">" {print c; c=0;printf substr($$0,2,100) "\t"; } \
			$$0 !~ ">" {c+=length($$0);} END { print c; }' $(REF) > $@

##################################
#### 4. COLLINEARITY ANALYSIS ####
##################################

# Preparing input file for collinearity analysis. Run MCScanX on files manually
.PHONY : collinearity
collinearity : $(RNA)/assembled/ $(CORRESP)
	rm -rf $(MCSX-IN)/
	mkdir -p $(MCSX-IN)/
	bash $(MCSX-SRC) -g $(RNA)/assembled/transcripts.gtf \
			-o $(MCSX-IN) \
			-r $(REF) \
			-c $(CORRESP) \
			$${LOCAL:+-l}
	bash src/circos_conf/circos_input_gen.sh $(REF)


# Assembling transcripts and measuring coverage along genome
$(RNA)/assembled/ :
	bash $(RNA-SRC) -a $(BAM) \
			-r $(OLD-REF) \
			-o $@ \
			$${LOCAL:+-l}

# Contig correspondance file between old (unanchored)
# and new (anchored) assembly
$(CORRESP):
	bash src/convert_coord/corresp_contigs.sh -O $(OLD-REF) \
			-g $(LOG) \
			-N $(REF) \
			-c $(CORRESP) \
			$${LOCAL:+-l} \
			2> $(LOG)/corresp.log

#########################################
#### 5. ANALYSIS OF WILD WGS SAMPLES ####
#########################################

$(WGS)/mapped : $(WGS)/raw
	bash src/wgs_wild/bwa_wgs.sh -d $(WGS) \
			-r $(REF) \
			$${LOCAL:+-l}


$(WGS)/variant/hap.wild.matrix : $(WGS)/mapped
	bash src/wgs_wild/snps_wgs.sh -d $(WGS) \
			-r $(REF) \
			$${LOCAL:+-l}

.PHONY : wgs_qc
wgs_qc : $(WGS)/mapped $(WGS)/raw
	bash src/wgs_wild/qc_gen.sh -d $(WGS) \
			-r $(REF) \
			$${LOCAL:+-l}

.PHONY : wgs_wild
wgs_wild : $(CORRESP) $(SIZES) $(WGS)/variant/hap.wild.matrix
	# Compute PI diversity in sliding windows
	Rscript src/wgs_wild/compute_PI.R \
		--in $(WGS)/variant/hap.wild.matrix.txt \
		--out $(WGS)/stats/win_w100_t10_PI.tsv \
		--mode 'window' \
		--step_size 10 \
		--win_size 100
	# Compute PI diversity per site
	Rscript src/wgs_wild/compute_PI.R \
		--in $(WGS)/variant/hap.wild.matrix.txt \
		--out $(WGS)/stats/sites_PI.tsv \
		--mode 'site'

####################
#### MISC RULES ####
####################

# Automatically download data from the paper and places it in the correct folders to
# run the analysis
.PHONY download
download: 
	wget https://zenodo.org/record/1488603/files/20181115_csd_data.tar.gz
	tar xzvf 20181115_csd_data.tar.gz
	rm 20181115_csd_data.tar.gz
	mkdir -p $(WGS)/raw
	bash $(MISC)/sra_autodl.sh SRX5006708 SRX5006729 $(WGS)/raw
	bash $(MISC)/sra_autodl.sh SRX5003361 SRX5003929 $(PROC)

$(WGS)/raw/%.fastq:
	wgs_prefix="SRX50067{08..29}"
	mkdir -p $(WGS)/raw
# Transform basepair positions for SNPs from the association mapping into 
# centiMorgan.
.PHONY: bp2cm
bp2cm:
	# Compute coordinates of linkage map markers in new genome
	# Compute cM/BP ratio in each inter-marker interval
	bash src/convert_coord/bp2cm.sh -r $(OLD-REF) -n $(REF) -o $(LINKMAP)/bp2cm/bp2cm \
									-m $(LINKMAP)/linkage_markers.csv
	# Apply transformation to CSD dataset to get cM coordinates
	bash src/convert_coord/apply_bp2cm.sh $(ASSOC)/case_control/case_control_all.tsv \
										  $(LINKMAP)/bp2cm/bp2cm.tsv \
										  $(LINKMAP)/bp2cm/bp2cm_csd_snps.tsv


# Check genome completeness using BUSCO
.PHONY : busco
busco :
	bash src/misc/genome_completeness.sh -r $(REF) \
			-d data/ref_genome/busco/

