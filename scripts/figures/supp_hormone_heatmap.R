# Figure of hormone heatmaps.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")
source("scripts/canola_pseudomonas_R_code.R")

library(cowplot)
library(ggplotify)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

ET_genes <- read.table("At_GO/hormone_genes/ET.txt", header=F, stringsAsFactors = FALSE)$V1
JA_genes <- read.table("At_GO/hormone_genes/JA.txt", header=F, stringsAsFactors = FALSE)$V1
SA_genes <- read.table("At_GO/hormone_genes/SA.txt", header=F, stringsAsFactors = FALSE)$V1

# Read in log2fold ratios.
ratios_AT <- read.table("At_deseq2_outfiles/At_homolog_deseq2_log2fold.txt",
                        header=TRUE, 
                        sep="\t",
                        stringsAsFactors = FALSE,
                        quote="",
                        comment.char = "")

rownames(ratios_AT) <- ratios_AT$At_gene

# Get info columns and remove from original df.
ratios_AT_info <- ratios_AT[, c(1, 2)]
ratios_AT <- ratios_AT[, -c(1, 2)]

# Make all names unique.
ratios_AT_info$unique_name <- make.unique(ratios_AT_info$NAME, sep = ".")

ratios_AT[ratios_AT > 2] <- 2
ratios_AT[ratios_AT < -2] <- -2

# Give columns clearer names
colnames(ratios_AT) <- c("Day 1 Shoot", "Day 3 Shoot", "Day 5 Shoot", "Day 1 Root", "Day 3 Root", "Day 5 Root")

breaksList = seq(-2, 2, by = 0.1)

SA_heatmap <- pheatmap(ratios_AT[SA_genes, ],
                       clustering_distance_rows = "euclidean",
                       clustering_method = "average",
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       treeheight_row = 0,
                       gaps_col=c(3, 3, 3),
                       angle_col=45,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                       breaks = breaksList,
                       main="Salicylic Acid Pathway")

JA_heatmap <- pheatmap(ratios_AT[JA_genes, ],
                       clustering_distance_rows = "euclidean",
                       clustering_method = "average",
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       treeheight_row = 0,
                       gaps_col=c(3, 3, 3),
                       angle_col=45,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                       breaks = breaksList,
                       main="Jasmonic Acid Pathway")
                       
ET_heatmap <- pheatmap(ratios_AT[ET_genes, ],
                       clustering_distance_rows = "euclidean",
                       clustering_method = "average",
                       cluster_rows = TRUE,
                       cluster_cols = FALSE,
                       treeheight_row = 0,
                       gaps_col=c(3, 3, 3),
                       angle_col=45,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                       breaks = breaksList,
                       main="Ethylene Pathway")





pdf(file = "plots/main/Supp_hormone_heatmaps.pdf", width=12, height=6, onefile=FALSE)

plot_grid(ggplot() + theme_void(), as.grob(SA_heatmap), as.grob(JA_heatmap), as.grob(ET_heatmap),
          nrow=1, rel_widths = c(0.1, 1, 1, 1), labels=c('', 'A', 'B', 'C'))

dev.off()

