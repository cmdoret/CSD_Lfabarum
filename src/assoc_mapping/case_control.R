# In this script I use a naive test for allelic association
# to test for association with homozygosity in males. The test should be performed 
# only within a category of family (inferred from diploid male production).
# Cyril Matthey-Doret
# 18.08.2017

#======== LOAD DATA =========####
# Loading libraries for data processing and visualisation
packs <- c("dplyr","readr","ggplot2","gridExtra","viridis", "ggrepel")
packs <- suppressPackageStartupMessages(sapply(packs, library, character.only=T, quietly=T))
in_args <- commandArgs(trailingOnly = T)
signif <- 0.00001
# Input table with proportion of homozygous individuals
hom_path <-in_args[1]
sum_stat <- read_tsv(file=hom_path, col_names=T, col_types = "icidddddd")
tig_sizes <- read_tsv(in_args[3], col_names=F, col_types = cols())
colnames(tig_sizes) <- c("Chr","size")
# Output folder
out_folder <- in_args[2]

#======= PROCESS DATA =======####
# Cleaning input table: removing SNPS absent in all samples
sum_stat <- sum_stat[sum_stat$N.Males>0 & sum_stat$N.Females>0,]
sum_stat <- sum_stat[sum_stat$Prop.Hom<1,]

#==== COMPUTE STATISTICS ====####
# Output single file with corrected significant p-values with
# associatied SNP and category.
# Abbreviations: M: Males, F: Females, T: Male+Female, 
# o:homozygous, e:heterozygous, t:hom+het, E: Expected

get_fisher <- function(df){
  # Computes fisher exact test on one row of the dataframe.
  # alternative=greater -> test if males are more heterozygous than females
  # alternative=less -> test if males are more homozygous than females
  mat <- matrix(as.numeric(df[c("Fo","Mo","Fe","Me")]), ncol=2)
  f <- fisher.test(as.table(mat), alt="less")
  return(f$p.value)
}

# Generating dataframe containing input values for the fisher test
# for a SNP in each row: Fo, Fe, Mo, Me. Original number of samples 
# in each category is obtained by multiplying the prop. of homozygosity 
# from the table by the total number of samples, and rounding to 
# closest integer.
odds_list <- sum_stat %>% 
  rename(Ft = N.Females, Mt = N.Males, Tt = N.Samples) %>%
  mutate(Fo = Ft * Prop.Hom.F, Mo = Mt * Prop.Hom.M, 
         Fe = Ft * (1-Prop.Hom.F), Me = Mt * (1-Prop.Hom.M),
         To = Tt * Prop.Hom, Te = Tt * (1-Prop.Hom)) %>%
  select(-Prop.Hom, -Prop.Hom.F, -Prop.Hom.M) %>%
  mutate_at(funs(round(.,0)), .vars = c("Fo","Fe","Mo","Me"))

# Applying the fisher test to each row. Each row gets converted
# into the 2x2 contingency table:
# Fo Fe
# Mo Me
# And one sided fisher test is performed to test if the odds ratio:
# (Fo/Fe)/(Mo/Me)
# is significantly less than 1.
odds_list$fisher <- apply(odds_list, 1,  get_fisher)

# Correcting p-values for multiple testing
odds_list$fisher <- p.adjust(odds_list$fisher, method = "BH")
# Transform to qvalue -log10 for convenience
odds_list$fisher <- round(-log10(odds_list$fisher),4)

#========= VISUALISE ========####
# Only including contigs that have been placed into chromosomes
odds_chrom <- odds_list[grep("chr.*",odds_list$Chr),]
# Unmapped contigs
odds_cont <- odds_list[grep("tig.*",odds_list$Chr),]

# Computing raw effect size:
# Proportion of homozygous males is a good measure for CSD
odds_chrom$effect_str <- odds_chrom$Me/(odds_chrom$Mo+odds_chrom$Me)
odds_chrom$effect_str[odds_chrom$effect_str>0.5] <- 0.5

# Saving manhattan plot to pdf
pdf(paste0(out_folder, "/../plots/","case_control_hits.pdf"), width=12, height=12)
ggplot(data=odds_chrom, aes(x=BP/1000000, y=fisher)) + 
  geom_point(aes(col=effect_str), size = 2) + 
  facet_grid(~Chr, space='free_x', scales = 'free_x') +  
  geom_hline(aes(yintercept=-log10(signif)), lty=2, col='red') + xlab("Genomic position (Mb)") + 
  ylab("-log10 p-value") + ggtitle("Case-control association test for CSD") + 
  ylim(c(0,10)) + theme_bw() + guides(col=guide_colorbar(title="Het. males")) + 
  scale_color_viridis()

dev.off()

viz_cont <- tig_sizes %>%
  arrange(-size) %>%
  filter(Chr %in% odds_cont$Chr) %>%
  mutate(cumSize = lag(cumsum(size), default=0), 
         idx = 1:n(), 
         colVal = idx %% 2 + 1) %>%
  inner_join(., odds_cont, by="Chr") %>%
  arrange(-size) %>%
  group_by(Chr) %>%
  mutate(sig_chr = ifelse(any(fisher > -log10(signif)), yes = 1, no = 0),
         show_name = case_when(fisher > -log10(signif) ~ Chr, TRUE ~ ""))

rects <- viz_cont %>% 
  filter(!duplicated(Chr), idx %% 2) %>% 
  mutate(start = cumSize, end = cumSize + size) %>%
  select(Chr, start, end)
zoom <- 10^6
ggplot(data=viz_cont, aes(x=(cumSize+BP)/zoom, y=fisher, col=as.factor(sig_chr))) +  
  geom_rect(data=rects, inherit.aes = F, aes(xmin=start/zoom, xmax=end/zoom, ymin=0, ymax=10), fill="grey80", alpha=0.4) +
  geom_point() +
  geom_hline(aes(yintercept=-log10(signif)), lty=2, col='red') + 
  geom_text_repel(aes(label=show_name),point.padding=2) +
  xlab("Genomic position [Mb]") + ylab("-log10 p-value") + 
  ggtitle("Case-control association test for CSD: Unordered contigs") + 
  ylim(c(0,10)) + theme_bw() + guides(col=F)

#======= WRITE OUTPUT =======####
# Writing table of significant hits to text file
odds_chrom_sig <- odds_chrom[odds_chrom$fisher>=-log10(signif),]
odds_cont_sig <- odds_cont[odds_cont$fisher>=-log10(signif),]
write.table(odds_chrom_sig, paste0(out_folder, "case_control_hits.tsv"), 
            sep='\t', row.names=F, quote=F)

write.table(odds_chrom, paste0(out_folder, "case_control_all.tsv"), 
            sep='\t', row.names=F, quote=F)

write.table(odds_cont_sig, paste0(out_folder, "unanchored_hits.tsv"), 
            sep='\t', row.names=F, quote=F)

write.table(odds_cont, paste0(out_folder, "unanchored_all.tsv"), 
            sep='\t', row.names=F, quote=F)

