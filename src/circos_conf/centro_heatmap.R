# This script fits a smooth curve to model the proportion of homozygosity along 
# chromosomes. This smoothing is performed with local regression (loess) and 
# values of the fitted curve are extracted at regular intervals (defined by resolution)
# These values are written into a circos-formatted table to produce a heatmap. 
# Cyril Matthey-Doret
# 24.11.2017

#==== SELECT parameters ====#
span_val <- 0.60  # Proportion of SNPs to be included in each local regression
map_res <- 10000  # Bin size for the heatmap
hom_path <- "data/assoc_mapping/grouped_outpool_prophom.tsv"
out_path <- "data/circos/centro.lf.txt"

#==== LOAD PACKAGES AND DATA ====#
# Loading packages for: Nice plots, data manipulation, colors, fast loading
pack <- c("ggplot2","dplyr","viridis", "readr")
lapply(pack, require, character.only = TRUE)
# Loading chromosomes sizes
chr_sizes <- read.table("data/circos/karyotype.lf.txt")
# Loading homozygosity table for RADseq SNPs
fix0 <- read_tsv(file=hom_path, col_names=T, col_types = "icidddddd")


#==== CLEAN AND PROCESS DATA ====#
# Formatting chromosome size table
chr_sizes <- chr_sizes[,c(7,6)]
colnames(chr_sizes) <- c("Chr","len")
# Filtering SNPs: present in at least 1 sample
fix <- fix0[fix0$N.Samples>0,]
# Filtering SNPs: Anchored in chromosome
fix <- fix[grep("chr.*",fix$Chr),]
# Creating empty output table for heatmap values
chr_sizes <- chr_sizes %>% 
  group_by(Chr, len) %>%
  summarise(N=ceiling(len/map_res))
# Heatmap data should be in the format: chr start end val
#heat <- data.frame(matrix(nrow=sum(chr_sizes$N), ncol=4))
#colnames(heat) <- c("Chr", "start", "end", "hom")
# Chromosome name for each window
#heat$Chr <- unlist(by(chr_sizes[,c("Chr","N")], INDICES = chr_sizes$Chr, 
#                      function(x) rep(x[["Chr"]],x[["N"]])),use.names = F)
# Adding coordinate of each window
#heat$start <- unlist(by(chr_sizes$len, INDICES = chr_sizes$Chr, 
#                        function(x) seq(0,max(x),by=map_res)), use.names = F)
#heat$end <- heat$start+(map_res-1)

#==== COMPUTE LOCAL REGRESSIONS ====#
mod <- data.frame(Chr= fix$Chr, start = numeric(nrow(fix)), end = numeric(nrow(fix)), hom = numeric(nrow(fix)))
for(chrom in unique(fix$Chr)){
  tmp_mod <- loess(data=fix[fix$Chr==chrom,], degree=2,
                 formula=Prop.Hom~BP, weights=N.Samples, span=span_val, model=T)
  tmp_start <- lag(tmp_mod$x)+1
  tmp_start[1] <- 0
  mod[mod$Chr==chrom,] <- data.frame(rep(chrom), tmp_start, tmp_mod$x, tmp_mod$fitted)
  mod$end[length(tmp_start)] <- chr_sizes$len[chr_sizes$Chr==chrom]
}
mod$hom <- round(mod$hom, 5)
#==== VISUALIZE ====#
# Will only run if the user chose one single value for each parameter.
#ggplot(data=mod, aes(x=start, y=hom, col=log(hom))) + geom_point() + 
#  facet_grid(~Chr, scales="free_x", space="free_x") + theme_bw() + scale_color_viridis()

#==== WRITE OUPUT ====#

write.table(mod, file=out_path, sep=' ', row.names = F, col.names = F, quote = F)
