# This script is used to differentiate between telocentric and metacentric
# chromosomes and locate each centromere. This is done by using the density 
# of heterozygous sites in the mother that become homozygous in offspring 
# along the chromosome.
# 30.07.2017
# Cyril Matthey-Doret

#==== SELECT PARAMETERS ====#
wsize_range <- 30  # Size of the moving average window
sp_range <- 0.40  # Proportion of SNPs to be included in each local regression
# wsize_range <- seq(5, 40, 1)
# sp_range <- seq(0.15, 1, 0.01)

in_path <- commandArgs(trailingOnly = T)[1]
ind_list <- commandArgs(trailingOnly = T)[2]
out_path <- commandArgs(trailingOnly = T)[3]

#==== LOAD PACKAGES AND DATA ====#
packs <- c("ggplot2","dplyr","viridis","zoo", "readr")
packs <- sapply(packs, function(x) 
  suppressPackageStartupMessages(library(x, quietly=T, character.only=T)))

indv <- read.table(ind_list, header=T)
fix0 <- read.table(in_path, header=T, na.strings='NA', sep='\t')
fix <- fix0[fix0$N.Samples>0,]
fix <- fix[grep("chr.*",fix$Chr),]
fix$Chr <- droplevels(fix$Chr)

#==== COMPUTE SLIDING MEANS ====#
# Allows to try different values of window size in a range
# Note that the step size is fixed to 1 but can be changed 
# directly in the code if neded.
pdf(paste0(out_path, "/plots/param_centro.pdf"))
virilist <- viridis(n=length(wsize_range))
colindex <- 1
par(mfrow=c(2,6))
centrolist <- list()
centrolist[['slideMean']] <- data.frame(pos = rep(0,6), Chr=levels(fix$Chr), val=rep(0,6))
store_means <- data.frame(BP=numeric(0),Chr=numeric(0),S.Mean=numeric(0))
for(chrom in levels(fix$Chr)){
  plot(x=c(),y=c(),xlim=c(0,max(fix$BP[fix$Chr==chrom])), ylim=c(-0.05,1),
       xlab=chrom,ylab="Homozygosity")
  abline(h=-0.01)
  for(w in wsize_range){
    HOM = zoo(fix$Prop.Hom[fix$Chr==chrom],order.by = fix$BP[fix$Chr==chrom])
    sliMean = rollapply(HOM, width=w, by=1, FUN=mean, partial=F)
    bp_idx = index(sliMean)
    points(bp_idx[order(bp_idx)],sliMean[order(bp_idx)], type='l',col=alpha(virilist[colindex],0.4),lwd=1.3)
    points(bp_idx[which(sliMean==min(sliMean, na.rm=T))],rep(-0.01,length(which(sliMean==min(sliMean, na.rm=T)))),
           col=alpha(virilist[colindex],0.4),pch=16, cex=1.5)
    centrolist$slideMean$pos[centrolist$slideMean$Chr==chrom] <- bp_idx[which(sliMean==min(sliMean, na.rm=T))]
    centrolist$slideMean$val[centrolist$slideMean$Chr==chrom] <- sliMean[sliMean==min(sliMean, na.rm=T)]
    colindex <- colindex+1
  }
  colindex <- 1
  tmp_df <- data.frame(BP=bp_idx, Chr=rep(chrom), S.Mean=sliMean)
  store_means <- rbind(store_means,tmp_df)
}

#==== COMPUTE LOCAL REGRESSIONS ====#

# Allows to try all different values of span in a range.
virilist <- viridis(n=length(sp_range))
colindex <- 1
#par(mfrow=c(1,6))
centrolist[['loess']] <- data.frame(pos = rep(0,6), Chr=levels(fix$Chr), val=rep(0,6))
for(chrom in levels(fix$Chr)){
  plot(x=c(),y=c(),xlim=c(0,max(fix$BP[fix$Chr==chrom])), ylim=c(-0.05,1),
       xlab=chrom,ylab="Homozygosity")
  abline(h=-0.01)
  for(sp in sp_range){
    mod <- loess(data=fix[fix$Chr==chrom,], degree=2,
                 formula=Prop.Hom~BP, weights=N.Samples, span=sp, model=T)
    points(mod$x[order(mod$x)],mod$fitted[order(mod$x)], type='l',col=alpha(virilist[colindex],0.4),lwd=1.3)
    points(mod$x[mod$fitted==min(mod$fitted)],rep(-0.01,length(mod$x[mod$fitted==min(mod$fitted)])),
           col=alpha(virilist[colindex],0.4),pch=16, cex=1.5)
    centrolist$loess$pos[centrolist$loess$Chr==chrom] <- mod$x[mod$fitted==min(mod$fitted)]
    centrolist$loess$val[centrolist$loess$Chr==chrom] <- min(mod$fitted)
    colindex <- colindex+1
  }
  colindex <- 1
}
dev.off()
#==== VISUALIZE ====#
# Will only run if the user chose one single value for each parameter.
if(length(sp_range)==1 & length(wsize_range)==1){
  zoomfactor <- 1000000  # For aesthetics
  LoessPlot <- ggplot(fix, aes(x=BP/zoomfactor, y=Prop.Hom, weight=N.Samples)) + facet_grid(~Chr, scales = 'free_x',space='free_x') + 
    geom_point(col='grey66',aes(size=N.Samples), alpha=0.4) + geom_smooth(aes(col='Local regression'), 
                                           method='loess', span=sp_range, fill='#EE9999') + 
    xlab("Genomic position (Mb)") + ylab("Homozygosity") + 
    geom_segment(data=centrolist$loess, aes(x=pos/zoomfactor,xend=pos/zoomfactor, y=0, yend=val, col='Local regression'), 
                 lty=2, lwd=1.1, inherit.aes = F) + theme_bw() +
    geom_line(data=store_means, aes(x=BP/zoomfactor, y=S.Mean, col='Moving average'), lwd=1.1, inherit.aes = F) + 
    geom_segment(data=centrolist$slideMean, aes(x=pos/zoomfactor,xend=pos/zoomfactor, y=0, yend=val, col='Moving average'), 
                 lty=2, lwd=0.9, inherit.aes = F) + scale_size_continuous(range = c(0,5))+
    scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") + 
    guides(color=guide_legend(title="Method",override.aes=list(fill=NA))) + 
    geom_point(data=centrolist$slideMean, aes(x=pos/zoomfactor, y=0, col='Moving average'), inherit.aes = F, size=2) + 
    geom_point(data=centrolist$loess, aes(x=pos/zoomfactor, y=0, col='Local regression'), inherit.aes = F, size=2)
  ggsave(filename = paste0(out_path, "/plots/final_centro.pdf"), plot = LoessPlot, height = 10, width = 14)
  write.table(centrolist$slideMean, file=paste0(out_path, "/centrolist.tsv"), sep='\t',row.names = F, quote = F)
}

