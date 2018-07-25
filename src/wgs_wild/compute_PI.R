#!/usr/bin/env Rscript
# Computing PI nucleotidic diversity along genome, both using a sliding window along the genome or per site at all 
# RADseq SNPs from the RADseq analysis.
# Cyril Matthey-Doret
# 23.03.2018

# Load libraries ####
packs <- c("dplyr","readr","Rcpp","optparse", "RcppRoll")
packs <- sapply(packs, function(x) suppressPackageStartupMessages(library(x, quietly=T, character.only=T)))

# Parse CL arguments ####
option_list = list(
  make_option(c("-i", "--in"), 
              type="character", 
              default=NA, 
              help="Input SNP matrix.", 
              metavar="character"),
  make_option(c("-o", "--out"), 
              type="character", 
              default=NA, 
              help="Output file where to store PI values.", 
              metavar="character"),
  make_option(c("-m", "--mode"), 
              type="character", 
              default="site", 
              help="Mode: Compute PI by 'window' or by 'site' [default %default]."),
  make_option(c("-s","--sites"), 
              type="character", 
              default=NA, 
              help="File containing a list of tab-separated chromosomes and positions (in base pairs) at which PI should be computed."),
  make_option(c("-w","--win_size"), 
              type="integer", 
              default=10000, 
              help="Size of windows in which PI is computed, in base pairs [default %default]."), 
  make_option(c("-t","--step_size"), 
              type="integer", 
              default=100, 
              help="Step between windows in which PI is computed, in base pairs [default %default].")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load data and handle arguments ####

# check output path is provided
if (!is.na(opt$out)) {
  # Loading genotype matrix. One column per haplotye, one row per SNP. Missing calls (.) as NA.
  out_path <- opt$out
} else stop("Output path not provided. See script usage (--help)")

# check "mode" is among accepted values
# check output path is provided
modes <- c("site","window")
if (!opt$mode %in% modes)
  stop(paste0("mode must be one of '", modes,"'. See script usage (--help)", collapse = ", "))

# check input file is provided
if (!is.na(opt$i)) {
  # Loading genotype matrix. One column per haplotye, one row per SNP. Missing calls (.) as NA.
  snp_tbl <- read_tsv(opt$i, col_names = F, na = ".", progress=F, col_types = cols())
} else stop("Input SNP matrix not provided. See script usage (--help)")

# Filter sites if not in window mode
if (!is.na(opt$sites)){
  if(opt$mode=='window') print("Ignoring site filtering: Not enabled in window mode.")
  else {
    # Filtering sites using the file provided via --sites
    site_tbl <- read_tsv(opt$s, col_names = F, progress = F, col_types = cols())
    snp_tbl <- snp_tbl %>% inner_join(., site_tbl,site_tbl,by=c("X1","X2"))
  }
}

# Getting genotypes into matrix format (without chr and pos)
snp_mat <- as.matrix(snp_tbl[,-c(1,2)], nrow=nrow(snp_tbl))

print("data loaded")

# Functions ####


# OBSOLETE (slow) Computes PI in a input genomic region. 
# Input: Each base pair is a row, each haplotype is a col.
cppFunction('double PIcpp_region(const NumericMatrix& seq){
            double mismatch = 0.0;
            int pairs = 0;
            int tot_alleles, n_other, n_alleles;
            NumericVector geno(seq.ncol()), freqs(max(na_omit(seq))+1);
            
            for(int s=0; s < seq.nrow(); s++){
            geno = seq( s, _);
            n_alleles = (max(na_omit(geno)) + 1);
            for(int allele = 0; allele < n_alleles; allele++){
            freqs[allele] = std::count(geno.begin(), geno.end(), allele);
            }
            tot_alleles = 0;
            for(int allele = 0; allele < n_alleles; allele++){
            tot_alleles += freqs[allele];
            }
            for(int allele = 0; allele < n_alleles; allele++){
            n_other = tot_alleles - freqs[allele];
            mismatch += (freqs[allele] * n_other);
            }
            pairs = pairs + tot_alleles * (tot_alleles - 1);
            }
            double pi = mismatch / pairs;
            return(pi);}'
)

# Pre compute mismatches and N pairs 
# Avoid redundant computations in overlapping windows later on
cppFunction('NumericVector prep_PI(const IntegerVector& geno){
  NumericVector v(2);
  if(all(is_na(geno))){
    v[0] = 0;
    v[1] = 0;
    return(v);}
  int n_alleles = (max(na_omit(geno)) + 1);
  int n_other;
  IntegerVector freqs(n_alleles);
  for(int allele = 0; allele < n_alleles; allele++){
    freqs[allele] = std::count(geno.begin(), geno.end(), allele);
  }
  int tot_alleles = sum(freqs);
  double mismatch = 0.0;
  for(int allele = 0; allele < n_alleles; allele++){
    n_other = tot_alleles - freqs[allele];
    mismatch = mismatch + freqs[allele] * n_other;
  }
  int pairs = tot_alleles * (tot_alleles - 1);
  v[0] = mismatch;
  v[1] = pairs;
  return(v);}'
)


# Computes PI in windows along the genome
slide_PI <- function(mat, pos, win, step){
  
  # Check whether parameters make sense
  if(any(win > max(pos), step == 0, win == 0, step > win)){
    print("Bad input parameters. Is your step size > window size, or matrix smaller than window ?")
    return()
  }
  
  # Computing only once N mismatches andx N pairs at all positions
  mat_stat <- matrix(ncol=2, byrow=T, apply(mat, FUN=prep_PI, MARGIN = 1))
  win_start <- seq(1, (max(pos) - (win-1)), step)
  n_win <- ceiling((max(pos)-(win-1))/step)
  # Initialize output structure to contain original index and pi
  out_mat <- matrix(nrow=n_win, ncol=3, byrow = F)
  out_mat[,1] <- win_start
  out_mat[,2] <- out_mat[,1] + win
  print(sprintf("Computing %i windows...", n_win))
  # Compute PI in all SNPs within window boundaries store pi
  mat_PI <- RcppRoll::roll_sum(mat_stat, n = win, by = step, align = 'left', na.rm = T)
  mat_PI <- mat_PI[seq(1,nrow(mat_PI),by=step),]
  out_mat[,3] <- mat_PI[,1]/mat_PI[,2]
  return(out_mat)
}


# Compute PI ####
out_PI <- snp_tbl[,1:2]
out_PI[,3] <- NA

out <- data.frame()
for(chr in unique(pull(snp_tbl, 1))){
  print(paste0("Computing PI for ", chr))
  chr_mat <- snp_mat[snp_tbl$X1==chr,]
  if(opt$m=="window"){
    # If positions are missing, add them as missing calls (easier for windows)
    full_chr <- matrix(nrow=max(snp_tbl[snp_tbl$X1==chr,2]),ncol=ncol(snp_mat))
    full_chr[snp_tbl$X2[snp_tbl$X1==chr],] <- chr_mat
    chr_slide <- slide_PI(full_chr, pos=snp_tbl$X2[snp_tbl$X1==chr], win=opt$win_size, step=opt$step_size)
    out <- rbind(out,cbind(rep(chr),chr_slide))
  } else{
    tmp_mat <- matrix(ncol=2, byrow=T, apply(chr_mat, FUN=prep_PI, MARGIN = 1))
    out_PI[snp_tbl$X1 == chr, 3] <- tmp_mat[,1]/tmp_mat[,2]
  }
}

print(paste0("Saving output to ", opt$o))
if(opt$m=="window") {out_PI <- out}
# Write PI values to output ####
write_tsv(out_PI, opt$out,col_names = F)
