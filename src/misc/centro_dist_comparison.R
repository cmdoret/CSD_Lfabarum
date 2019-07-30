# Comparing distance from centromeres between candidate SNPs and non-candidate SNPs
library('tidyverse')
library('gridExtra')
library('ggpubr')
snps <- read_delim('data/assoc_mapping/case_control/case_control_all.tsv', delim='\t')
centro <- read_delim('data/assoc_mapping/centro/centrolist.tsv', delim='\t')
head(centro)

snps_centro <- snps %>%
    inner_join(centro, by="Chr") %>%
    filter(Mt > 0 & Ft > 0) %>%
    mutate(centrodist = abs(BP - pos),
           candidate = ifelse(fisher >= 5, T, F)) %>%
    select(Chr, centrodist,Tt, fisher, candidate, effect_str)

zoomfactor <- 1000000

p1 <- ggplot(snps_centro, aes(x=centrodist / zoomfactor)) + 
    geom_histogram(data=snps_centro %>% filter(candidate == F)) + 
    geom_vline(data=snps_centro %>% filter(candidate == T), aes(xintercept=centrodist / zoomfactor), col='red')  +
    facet_grid(~Chr, space='free_x', scales='free_x') + 
    theme_classic() + 
    xlab("Absolute distance from centromere [Mbp]") +
    ylab("Number of SNPs")

p2 <- ggplot(snps_centro, aes(x=as.character(candidate), y=centrodist / zoomfactor)) + 
    geom_boxplot() +
    geom_point() + 
    stat_compare_means(comparisons = list(c("FALSE", "TRUE")), label = "p.signif")+ # Add significance levels
    stat_compare_means(label.x=0.5, label.y = 8*10e5 / zoomfactor) +
    theme_bw() +
    coord_flip() +
    xlab("CSD candidate") +
    ylab("Absolute distance to centromere [Mbp]")

pdf('data/assoc_mapping/centro/plots/centro_dist.pdf', width=16, height=6)
grid.arrange(nrow=2, p1, p2)
dev.off()