# Figure of concordance between qPCR and RNA-seq. Can make two scatterplots: one for root and one for shoot.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

library(cowplot)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Read in qPCR results and RNA-seq (lfcshrink version)
# ABCG40 = BnaA06g10230D and BnaC05g11890D
# BBE4 = BnaC05g20370D and BnaA09g28890D
# CYP710A1 = BnaC04g10740D
# RBCSF1 = All genes that match AT5G38420 (which is called "2B" in A. thaliana rather than "F1", but was the top BLAST hit)

# Define Bnapus gene sets to use.
gene_to_ids <- list()
gene_to_ids[["ABCG40"]] <- c("BnaA06g10230D", "BnaC05g11890D")
gene_to_ids[["BBE4"]] <- c("BnaC05g20370D", "BnaA09g28890D")
gene_to_ids[["CYP710A1"]] <- "BnaC04g10740D"

Bnapus_deseq2_out <- read.table("Bnapus_deseq2_outfiles/deseq2_log2fold.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", row.names=1)
gene_to_ids[["RBCSF1"]] <- rownames(Bnapus_deseq2_out)[which(Bnapus_deseq2_out$Athaliana_top_hit == "AT5G38420")]

qPCR <- read.table("tables/qPCR_data.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
qPCR$day_str <- paste("day", as.character(qPCR$day), sep="")

qPCR$root_qPCR_se <- qPCR$root_qPCR_sd / sqrt(3)
qPCR$shoot_qPCR_se <- qPCR$shoot_qPCR_sd / sqrt(3)

rna_levels <- list()

for(tissue in c("root", "shoot")) {
  rna_levels[[tissue]] <- list()
  
  for(day in c("day1", "day3", "day5")) {
    
    rna_levels[[tissue]][[day]] <- data.frame(matrix(NA, nrow=4, ncol=6))
    colnames(rna_levels[[tissue]][[day]]) <- c("Gene", "Day", "rnaseq_mean", "qpcr_mean", "rnaseq_se", "qpcr_se")
    rownames(rna_levels[[tissue]][[day]]) <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
    rna_levels[[tissue]][[day]]$Gene <- c("ABCG40", "BBE4", "CYP710A1", "RBCSF1")
    rna_levels[[tissue]][[day]]$Day <- day
    
    if(tissue == "root") {
      lfcshrink_file <- paste("Bnapus_deseq2_outfiles/", day, "_results_RC_RI_shrink.txt", sep="")
      qPCR_mean_col <- "root_qPCR_mean"
      qPCR_se_col <- "root_qPCR_se"
    } else if(tissue == "shoot") {
      lfcshrink_file <- paste("Bnapus_deseq2_outfiles/", day, "_results_SC_SI_shrink.txt", sep="")
      qPCR_mean_col <- "shoot_qPCR_mean"
      qPCR_se_col <- "shoot_qPCR_se"
    }
    
    tmp_lfcshrink <- read.table(lfcshrink_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
    
    for(gene in names(gene_to_ids)) {
      rna_levels[[tissue]][[day]][gene, "rnaseq_mean"] <- mean(tmp_lfcshrink[gene_to_ids[[gene]], "log2FoldChange"])
      rna_levels[[tissue]][[day]][gene, "rnaseq_se"] <- mean(tmp_lfcshrink[gene_to_ids[[gene]], "lfcSE"])
    }
    
    rna_levels[[tissue]][[day]]$qpcr_mean <-qPCR[which(qPCR$day_str == day), qPCR_mean_col]
    rna_levels[[tissue]][[day]]$qpcr_se <- qPCR[which(qPCR$day_str == day), qPCR_se_col]
  }
}

root_rna <- do.call("rbind", rna_levels$root)
shoot_rna <- do.call("rbind", rna_levels$shoot)

root_rna[which(root_rna$Day == "day1"), "Day"] <- "Day 1"
root_rna[which(root_rna$Day == "day3"), "Day"] <- "Day 3"
root_rna[which(root_rna$Day == "day5"), "Day"] <- "Day 5"

root_scatterplot <- ggplot(root_rna, aes(x=rnaseq_mean, y=qpcr_mean, fill=Day, shape=Gene)) + 
                            geom_point(size=5, color="black") +
                            scale_fill_manual(values=c("white", "black", "grey")) +
                            scale_shape_manual(values=c(21, 22, 23, 24)) +
                            guides(fill = guide_legend(override.aes=list(shape=22))) +
                            geom_errorbar(aes(ymin=qpcr_mean - qpcr_se,
                                              ymax=qpcr_mean + qpcr_se),
                                          width=0.02,
                                          size=0.5) +
                            geom_errorbarh(aes(xmin=rnaseq_mean - rnaseq_se,
                                               xmax=rnaseq_mean + rnaseq_se),
                                           size=0.5) +
                            geom_vline(xintercept = 0, color="grey", linetype="dashed", size=0.5) +
                            geom_hline(yintercept = 0, color="grey", linetype="dashed", size=0.5) +
                            ggtitle("Root") +
                            ylab(expression('RT-qPCR log'[2]*'-fold change')) +
                            xlab(expression('RNA-seq log'[2]*'-fold change')) +
                            xlim(-5, 10) +
                            ylim(-5, 10) +
                            theme(legend.justification = c(0.05, 1), legend.position = c(0.05, 1),
                                  legend.background = element_blank(),
                                  legend.box.background = element_rect(colour = "black", fill="white"),
                                  plot.title = element_text(hjust = 0.5))


shoot_rna[which(shoot_rna$Day == "day1"), "Day"] <- "Day 1"
shoot_rna[which(shoot_rna$Day == "day3"), "Day"] <- "Day 3"
shoot_rna[which(shoot_rna$Day == "day5"), "Day"] <- "Day 5"

shoot_scatterplot <- ggplot(shoot_rna, aes(x=rnaseq_mean, y=qpcr_mean, fill=Day, shape=Gene)) + 
  geom_point(size=5, color="black") +
  scale_fill_manual(values=c("white", "black", "grey")) +
  scale_shape_manual(values=c(21, 22, 23, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=22))) +
  geom_errorbar(aes(ymin=qpcr_mean - qpcr_se,
                    ymax=qpcr_mean + qpcr_se),
                width=0.02,
                size=0.5) +
  geom_errorbarh(aes(xmin=rnaseq_mean - rnaseq_se,
                     xmax=rnaseq_mean + rnaseq_se),
                 size=0.5) +
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=0.5) +
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=0.5) +
  ggtitle("Shoot") +
  ylab(expression('RT-qPCR log'[2]*'-fold change')) +
  xlab(expression('RNA-seq log'[2]*'-fold change')) +
  xlim(-5, 10) +
  ylim(-5, 10) +
  theme(legend.justification = c(0.05, 1), legend.position = c(0.05, 1),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", fill="white"),
        plot.title = element_text(hjust = 0.5))


pdf(file = "plots/main/Figure4_qPCR_vs_RNAseq.pdf", width=12, height=6, onefile=FALSE)
plot_grid(root_scatterplot, shoot_scatterplot,
          nrow=1,
          labels=c('A', 'B'))
dev.off()

tiff(file = "plots/main/Figure4_qPCR_vs_RNAseq.tiff", width=12, height=6,
     compression = "lzw", res=600, units="in")
plot_grid(root_scatterplot, shoot_scatterplot,
          nrow=1,
          labels=c('A', 'B'))
dev.off()