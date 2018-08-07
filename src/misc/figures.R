### LOAD FILES
library(tidyverse); library(ggrepel); library(viridis); library(hexbin); library(zoo)

setwd('data/paper_output')
ploid <- read_tsv('d20_fixed.tsv', col_names=T)
indv <- read_tsv('individuals.tsv', col_names=T)
fix0 <- read_tsv('grouped_outpool_prophom.tsv', col_names=T)
tig_sizes <- read_tsv('tig.sizes.txt', col_names=F)
assoc_snp <- read_tsv('assoc_mapping/case_control_all.tsv', col_names=T)
assoc_tig <- read_tsv('assoc_mapping/unanchored_all.tsv', col_names=T)
win_pi <- read_tsv('wgs/win_w100_t10_PI.tsv', col_names = F) %>%
    rename(Chr = X1, Start = X2, End = X3, Pi = X4)

signif <- 10^-5
zoom <- 10^6



### Fig. 2: Manhattan plot
# Note: fisher variable = log10 of BH-corrected p-values from the fisher exact test

pdf("pdf/fig2_association.pdf")
ggplot(data=assoc_snp, aes(x=BP/zoom, y=fisher)) + 
  geom_point(size = 2) + 
  facet_grid(~Chr, space='free_x', scales = 'free_x') +  
  geom_hline(aes(yintercept=-log10(signif)), lty=2, col='red') + xlab("Genomic position (Mb)") + 
  ylab("-log10 p-value") + ggtitle("Case-control association test for CSD") + 
  ylim(c(0,10)) + theme_bw() +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
dev.off()

### Fig. S1: Ploidy separation
combined <- ploid %>% select(-Sex)
pdf("pdf/figS1_ploidy.pdf")
ggplot(data=ploid, aes(x=HOM, fill=Sex)) +
    geom_histogram(data=combined, fill='grey', binwidth = 0.01) +
    geom_histogram(binwidth = 0.01) +
    scale_fill_brewer(palette="Set1") +
    geom_vline(aes(xintercept=0.9), col='red', lty=2) +
    theme_minimal() +
    facet_grid(Sex~.) +
    xlab("Proportion of homozygous SNPs") +
    ylab("Number of individuals")
dev.off()

### Fig. S2: Hits on unanchored contigs

viz_cont <- tig_sizes %>%
  rename(Chr=X1, size=X2) %>%
  arrange(-size) %>%
  filter(Chr %in% unique(assoc_tig$Chr)) %>%
  mutate(cumSize = lag(cumsum(size), default=0), 
         idx = 1:n(), 
         colVal = idx %% 2 + 1) %>%
  inner_join(., assoc_tig, by="Chr") %>%
  arrange(-size) %>%
  group_by(Chr) %>%
  mutate(sig_chr = ifelse(any(fisher > -log10(signif)), yes = 1, no = 0),
         show_name = case_when(fisher > -log10(signif) ~ Chr, TRUE ~ ""))

rects <- viz_cont %>%
  filter(!duplicated(Chr), idx %% 2) %>% 
  mutate(start = cumSize, end = cumSize + size) %>%
  select(Chr, start, end)

pdf("pdf/figS2_unanchored.pdf")
ggplot(data=viz_cont, aes(x=(cumSize+BP)/zoom, y=fisher, col=as.factor(sig_chr))) +  
  geom_rect(data=rects, inherit.aes = F, aes(xmin=start/zoom, xmax=end/zoom, ymin=0, ymax=10), fill="grey80", alpha=0.4) +
  geom_point() +
  geom_hline(aes(yintercept=-log10(signif)), lty=2, col='red') + 
  geom_text_repel(aes(label=show_name),nudge_x=sample(1:2),nudge_y=sample(0.3:0.5), point.padding = 0.5) +
  xlab("Genomic position [Mb]") + ylab("-log10 p-value") + 
  ggtitle("Case-control association test for CSD: Unordered contigs") + 
  ylim(c(0,10)) + theme_bw() + guides(col=F) +
  scale_color_brewer(palette="Set1")
dev.off()

### Fig. S3: Manhattan with male heterozygosity and centimorgans

pdf("pdf/figS3_cM.pdf")
ggplot(data=assoc_snp, aes(x=cM, y=fisher)) + 
  geom_point(aes(col=effect_str), size = 2) + 
  facet_grid(~Chr, space='free_x', scales = 'free_x') +  
  geom_hline(aes(yintercept=-log10(signif)), lty=2, col='red') + xlab("Genetic distance (cM)") + 
  ylab("-log10 p-value") + ggtitle("Case-control association test for CSD") + 
  ylim(c(0,10)) +
  theme_bw() +
  guides(col=guide_colorbar(title="Het. males")) + 
  scale_color_viridis()
dev.off()

### Fig. S4: Centromeres location
wsize_range <- 50
sp_range <- 0.40

fix <- fix0[fix0$N.Samples>0,]
fix$Chr <- as.factor(fix$Chr)
fix <- fix[grep("chr.*",fix$Chr),]
fix$Chr <- droplevels(fix$Chr)

# Sliding means
centrolist <- list()
centrolist[['slideMean']] <- data.frame(pos = rep(0,6), Chr=levels(fix$Chr), val=rep(0,6))
store_means <- data.frame(BP=numeric(0),Chr=numeric(0),S.Mean=numeric(0))
for(chrom in levels(fix$Chr)){
  for(w in wsize_range){
    HOM = zoo(fix$Prop.Hom[fix$Chr==chrom],order.by = fix$BP[fix$Chr==chrom])
    sliMean = rollapply(HOM, width=w, by=1, FUN=mean, partial=F)
    bp_idx = index(sliMean)
    centrolist$slideMean$pos[centrolist$slideMean$Chr==chrom] <- bp_idx[which(sliMean==min(sliMean, na.rm=T))]
    centrolist$slideMean$val[centrolist$slideMean$Chr==chrom] <- sliMean[sliMean==min(sliMean, na.rm=T)]
  }
  tmp_df <- data.frame(BP=bp_idx, Chr=rep(chrom), S.Mean=sliMean)
  store_means <- rbind(store_means,tmp_df)
}

# Local regressions

# Allows to try all different values of span in a range.
#par(mfrow=c(1,6))
centrolist[['loess']] <- data.frame(pos = rep(0,6), Chr=levels(fix$Chr), val=rep(0,6))
for(chrom in levels(fix$Chr)){
  for(sp in sp_range){
    mod <- loess(data=fix[fix$Chr==chrom,], degree=2,
                 formula=Prop.Hom~BP, weights=N.Samples, span=sp, model=T)
    centrolist$loess$pos[centrolist$loess$Chr==chrom] <- mod$x[mod$fitted==min(mod$fitted)]
    centrolist$loess$val[centrolist$loess$Chr==chrom] <- min(mod$fitted)
  }
}

#==== VISUALIZE ====#
# Will only run if the user chose one single value for each parameter.
zoomfactor <- 1000000  # For aesthetics
LoessPlot <- ggplot(fix, aes(x=BP/zoomfactor, y=Prop.Hom, weight=N.Samples)) +
    facet_grid(~Chr, scales = 'free_x',space='free_x') + 
    geom_point(col='grey66',aes(size=N.Samples), alpha=0.4) +
    geom_smooth(aes(col='Local regression'), 
    method='loess', span=sp_range, fill='#EE9999') + 
    xlab("Genomic position (Mb)") + ylab("Homozygosity") + 
    geom_segment(data=centrolist$loess, 
                 aes(x=pos/zoomfactor,xend=pos/zoomfactor, 
                     y=0, yend=val, col='Local regression'), 
                 lty=2, lwd=1.1, inherit.aes = F) + 
    theme_bw() +
    geom_line(data=store_means, aes(x=BP/zoomfactor, y=S.Mean, col='Moving average'), lwd=1.1, inherit.aes = F) + 
    geom_segment(data=centrolist$slideMean, aes(x=pos/zoomfactor,xend=pos/zoomfactor, 
                                                y=0, yend=val, col='Moving average'), 
                 lty=2, lwd=0.9, inherit.aes = F) + scale_size_continuous(range = c(0,5))+
    scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") + 
    guides(color=guide_legend(title="Method",override.aes=list(fill=NA))) + 
    geom_point(data=centrolist$slideMean, aes(x=pos/zoomfactor, y=0, col='Moving average'), inherit.aes = F, size=2) + 
    geom_point(data=centrolist$loess, aes(x=pos/zoomfactor, y=0, col='Local regression'), inherit.aes = F, size=2)

pdf('pdf/figS4_centromeres.pdf')
LoessPlot
dev.off()

### Fig. S5: Recombination rate along genome
assoc_snp <- assoc_snp %>%
    group_by(Chr) %>%
    arrange(BP) %>%
    mutate(idx = 1:n())

pdf("pdf/figS5_recombination.pdf")
ggplot(data=assoc_snp, aes(x=idx, y=cM)) +
    geom_point() + 
    theme_bw() +
    facet_grid(~Chr) +
    xlab("SNP rank") + 
    ylab("Genetic (cM)")
dev.off()

### Fig. S6: PI nucleotidic diversity

hit_region <- assoc_snp %>%
    filter(fisher >= 5) %>%
    mutate(left  = BP - 500, 
           right = BP + 500)

pdf("pdf/figS6_pi.pdf")
ggplot(data=win_pi, aes(x=(Start + End)/(2*zoom), y=Pi)) +
    facet_wrap(~Chr, scales='free_x') +
    stat_bin_hex(bins=100) +
    theme_bw() +
    scale_fill_viridis() +
    geom_rect(data=hit_region, aes(xmin=left/zoom, xmax=right/zoom, 
                                   ymin=0, ymax=max(win_pi$Pi, na.rm=T)), 
              fill="red", alpha=0.3, inherit.aes = F) +
    xlab("Genomic position [Mb]") +
    ylab("Pi nucleotidic diversity (100bp windows)")
dev.off()

pval <- c()
hit_win <-  data.frame()
for(hit in 1:nrow(hit_region)){
    overlap <- win_pi[win_pi$Chr == pull(hit_region[hit, 'Chr']) &
                      abs(((win_pi$Start + win_pi$End) / 2) - pull(hit_region[hit, 'BP'])) < 500,]
    hit_win <- bind_rows(hit_win, overlap)

    pval <- append(pval, wilcox.test(overlap$Pi,
                      win_pi$Pi[win_pi$Chr == pull(hit_region[hit, 'Chr'])])$p.value)
}
hit_region$pval<- pval

ggplot(data=win_pi, aes(x=log10(Pi))) + 
    geom_histogram(binwidth=0.01, aes(y=..density..)) + 
    theme_bw()+
    geom_histogram(data=hit_win, binwidth=0.01, 
                   aes(x = log10(Pi), y = -..density..), fill='red') +
    facet_wrap(~Chr)
