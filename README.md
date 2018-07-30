
# Genetic basis of sex determination in a parasitoid wasp species

### Cyril Matthey-Doret, Casper Van der Kooi and Tanja Schwander

---
In this project, we use restriction-site associated DNA-sequencing (RAD-seq) and build a custom pipeline to locate and identify the complementary sex determination (CSD) locus/loci in the parasitoid wasp species _Lysiphlebus fabarum_. This repository contains a pipeline to map the reads using BWA and build loci using the different components of STACKS with optimal parameters. It was designed to run on a distributed High Performance Computing (HPC) environment with LSF but can also be ran on a personal machine.

## Instructions:

### Run the pipeline

To run the pipeline with the data provided:

1. Download or clone this repository.
2. Download the `data` archive (not available yet) into the main folder
3. ```cd``` to the main folder
4. Untar the data using ```tar -xzvf data.tar.gz data```
5. Run the  pipeline
    1. On a cluster with LSF:
        + `make` to run the STACKS pipeline
        + `make assoc_mapping` to run the association mapping (needs STACKS data)
        + `make collinearity` to compute collinearity blocks (needs STACKS data)
        + `make wgs_wild` to run the analysis of wild WGS samples
    2. On a local machine
        + `make ref_local` to run the STACKS pipeline
        + `make assoc mapping` to run the association mapping (needs STACKS data)
        + `make collinearity LOCAL=yes` to compute collinearity blocks (needs STACKS data)
        + `make wgs_wild LOCAL=yes` to run the analysis of wild WGS samples

### Use different input data

To run the STACKS pipeline with new data in the form of demultiplexed, trimmed single end reads in compressed fastq files (.fq.gz):

1. Describe your samples by writing 2 files named `popmap.tsv` and `individuals.tsv`, respectively. The structure of the `popmap.tsv` file is described on the [official STACKS website](http://catchenlab.life.illinois.edu/stacks/manual/) (here, populations should be the sex of individuals). The `individuals.tsv` file is a tab delimited text file with 4 columns with the following headers __included__:
    + Name: The names of samples. This should be the name of their data files (e.g. if the sample name is SAMPLE1, the corresponding reads file should be named SAMPLE1.fq.gz).
    + Sex: F for females and M for males.
    + Family: Clutches to which the individual belongs. These can be any combination of alphanumeric characters.
    + Generation: Useful if there are mothers and offspring. Values should be F3 for mothers and F4 for offspring.
2. Create an empty folder named data and place the 2 files inside. This folder needs to be located inside the same directory as src.
3. Place your (trimmed, demultiplexed) reads in a subfolder of data named `processed` and your reference genome in a subfolder named `ref_genome`. You will also need to edit the `REF` path in `config.mk` accordingly. If you wish to use different folder names, just edit the corresponding paths in `config.mk`.

4. Set the variable `D` in `config.mk` to 20 (minimum locus depth for STACKS populations). Type `make` in the command line (or `make ref_local` if running on a local machine). Once the pipeline has finished running, set the variable `D` back to 5 and type `make ploidy` to infer ploidy from the homozygosity of variant sites. Note the threshold selected to define ploidy is adapted to the dataset presented here. You might want to define a threshold yourself by inspecting the distribution of homozygosity (HOM variable) in `data/ploidy/thresholds/fixed.tsv` and the variable `HOM_PLOID` in `config.mk` to this value. Once you have modified the threshold, run `make ploidy` again to update the ploidy classification with the new threshold.
5. Type `make -B` to run the pipeline again without haploids.

### Dependencies:

The easiest and recommended way to install the correct versions of python and R packages is to install [Anaconda](https://www.anaconda.com/download) for python 2 and [setup an environment](https://conda.io/docs/user-guide/tasks/manage-environments.html). For convenience, a YAML environment file is provided in the `docs` folder and you can setup the conda environment directly from it using `conda env create -f docs/csd_env_conda.yml`. You can then activate it using `source activate csd_env` and all python and R dependencies should be enabled.

#### Core:

These dependencies are required to run the main analysis. This include the RADseq pipeline and association mapping.


* [FastQC 0.11.5](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): Quality control of sequencing data.
* [BWA 0.7.15](http://bio-bwa.sourceforge.net/)
* [STACKS 1.48](http://catchenlab.life.illinois.edu/stacks/) ([Download link](http://catchenlab.life.illinois.edu/stacks/source/stacks-1.48.tar.gz))
* [SAMtools 1.4](http://samtools.sourceforge.net/)
* [VCFtools 0.1.15](https://vcftools.github.io/)
* [BEDtools 2.26](http://bedtools.readthedocs.io/)
* [Trimmomatic 0.36](http://usadellab.org/cms/?page=trimmomatic)
* [R 3.4.3](https://www.r-project.org/)
  + [tidyverse 1.2.1](https://www.tidyverse.org)
  + [zoo 1.8.0](https://cran.r-project.org/web/packages/zoo/index.html)
  + [Rcpp 0.12.14](https://cran.r-project.org/web/packages/Rcpp/)
  + [RcppRoll 0.2.2](https://cran.r-project.org/web/packages/RcppRoll/)
  + [viridis 0.4.0](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
  + [argparse 1.1.0](https://cran.r-project.org/web/packages/argparse/index.html)
  + [gridExtra 2.3](https://cran.r-project.org/web/packages/gridExtra/index.html)
* [Python 2.7.x](https://www.python.org/)
  + [numpy 1.11](http://www.numpy.org/)
  + [pandas 0.19](http://pandas.pydata.org/)
  + [matplotlib 1.5](https://matplotlib.org/)
  + [biopython 1.72](http://biopython.org/)
  + [pybedtools 0.7.10](http://daler.github.io/pybedtools/)

#### Optional:

These tools are required only to run additional analyses.

* [DeepTools 2.4.2](http://deeptools.readthedocs.io/)
* [CuffLinks 2.2.1](http://cole-trapnell-lab.github.io/cufflinks/)
* [ncbi-BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [MCScanX](http://chibba.pgml.uga.edu/mcscan2/)
* [circos](http://circos.ca/)
* [BUSCO 3.0.2](https://busco.ezlab.org/)

### Scripts

The `src` folders contains all scripts required to run the analysis along with other programs used for automated report generation and benchmarking. Those scripts are organized into several sub folders:

* `assoc_mapping`: This folder contains scripts used to locate candidate CSD region(s).
  + `case_control.R`: Performs the actual association mapping, using the Fisher exact test.
  + `chrom_types.R`: Modeling the proportion of recombinant offspring along chromosomes to locate centromeres.
  + `process_genomic.py`: Process genomic output from STACKS' populations module by transforming numeric encoding of genotypes into homozygous/heterozygous/missing, removing loci that are either homozygous or missing in mothers from their families and computing proportion of homozygous individuals per sex/family at each site.

* `mapping`: This folder contains scripts used to map processed sequencing reads to the reference genome of _Lysiphlebus fabarum_.
  + `bwa_script.sh`: Coordinates the mapping of all samples using BWA, sending the output to `split_sam.pl`.
  + `split_sam.pl`: Parses the output sam files to split single hits and multiple hits into separate files and convert the files into bam format.

* `misc`: This folder contains various scripts used at different steps of the pipeline.
  + `parse_VCF.sh`: Uses vcftools to compute several statistics from the output VCF file returned by the populations module of STACKS and store them in text files inside the `vcftools` subfolder. This script is used in the `ploidy` Makefile rule, when inferring the ploidy of males.

* `ploidy`: This folder contains scripts required to classify males as diploid or haploid based on the genomic data.
  + `haplo_males.py`: Uses a threshold to infer the ploidy of males.
  + `blacklist_haploloci.py`: Generates a list of loci to blacklist in subsequent runs. Loci blacklisted are those heterozygous in more than 50% of haploid males, by default.

* `process_reads`: This folder contains the script used for demultiplexing, trimming and removing adaptors from raw sequencing reads. These are not implemented in the pipeline as the processed reads should be provided.

* `stacks_pipeline`: This folder contains scripts required to run the different components of the STACKS suite
  + `pstacks_script.sh`: Produces 'stacks' from processed reads, using the Pstacks module.
  + `cstacks_script.sh`: Constructs a catalogue of loci from the output of `pstacks_script.sh`. Only files containing at least 10% of the mean number of RAD-tags (computed over all files) are included in the catalogue to remove poor quality samples.
  + `group_sstacks.sh`: Copy pstacks and cstacks output files to the sample folder to provide a working directory for `multi_sstacks.sh`.
  + `sstacks_script.sh`: Produces 'match' files from pstacks and cstacks output files (stacks and catalogue, respectively).
  + `populations.sh`: Uses the populations module to compute populations statistics and generate different outputs from sstacks output files.


### Data files

Once the `data.tar.gz` has been uncompressed, the data folder should contain the following files:

* `annotations`: Contains genome annotations from the BIPAA.
* `individuals.tsv`: Detailed characteristic of each individuals: Name, Sex, Family and Generation where F4 are son/daughter and F3 is the mother.
* `ploidy`: contains information about the ploidy of individuals in the dataset. Ploidy information is stored in `thresholds/fixed.tsv.
* `popmap.tsv`: Population map required by STACKS to match sample names to population group (i.e. male and female).
* `processed`: This folder contains the processed RAD-tags generated for each sample using process_radtags.
* `ref_genome`: This folder contains the reference genome.
* `wgs`: contains the raw reads from whole genome sequencing of a wild population of _L. fabarum_ and `wgs_samples.tsv`, a list of sample names with their sex.

After the pipeline has been running, all intermediary and final output files will be generated and stored in their respective sub-folders inside `data`. Notable output files include:

* Association mapping hits: `data/assoc_mapping/case_control/case_control_hits.tsv`
* STACKS populations output files: `data/popula`
* Locations of centromeres: `data/assoc_mapping/centro`
* List of samples including ploidy and coverage statistics: `data/ploidy/thresholds/fixed.tsv`
* Nucleotidic diversity stats from WGS samples: `data/wgs_wild/stats/`
* List of SNPs with centiMorgan information extracted from the linkage map: `data/linkage_map/bp2cm/bp2cm_csd_snps.tsv`
