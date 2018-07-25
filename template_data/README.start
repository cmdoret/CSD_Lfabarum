# Reproducing the CSD analysis with different data

In the following instructions, the root folder is consedered to be the git repository. That is the folder which contains src.

0. run `git clone https://github.com/cmdoret/CSD_lfabarum.git/ and create a data subfolder `data` inside the newly created folder
1. Replace informations in `template_data/individuals.tsv` and `template_data/popmap.tsv` by actual values and move them into `data`
2. Move your fastq files to `data/processed`, respecting example naming scheme
3. Insert genome fasta file in `data/ref_genome`
5. Have a look at variables in config.mk and edit variables as necessary (e.g. path to reference genome, STACKS filtering parameters, ...). Make sure that D=20 (min STACKS populations depth for including loci is 20X).
6. run `make` if running on a cluster with LSF Otherwise, run `make LOCAL=TRUE`.
7. When the run has finished, run "make ploidy" to infer ploidy.
8. In config.mk, set D=5.
9. repeat step 6.
10. To run the association mapping, run `make assoc_mapping`
