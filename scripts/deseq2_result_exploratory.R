library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

setwd("/home/gavin/projects/pseudomonas/canola_pseudomonas_RNAseq/")

sample_deseq2 <- readRDS("At_deseq2/RDS_files/sample_deseq2.rds")

# Run variance stabilizing transformation blindly to compare samples for QA.
vst_blind <- vst(sample_deseq2, blind=TRUE)

# Plot sample distances.
sampleDists <- dist(t(assay(vst_blind)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vst_blind$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

qual_col <- c( "#a6cee3" ,  "#1f78b4" , "#b2df8a" , "#33a02c" , "#fb9a99" , "#e31a1c" ,
               "#fdbf6f" , "#ff7f00" , "#cab2d6" , "#6a3d9a" , "#ffff99" , "#b15928" )

# Plot PCA for same samples.
pcaData <- plotPCA(vst_blind, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) + 
  coord_fixed() +
  scale_color_manual(values=qual_col)

# Exploratory plots: (example)
day5_results_SC_SI_shrink <- readRDS("At_deseq2/RDS_files/day5_results_SC_SI_shrink.rds")
plotMA(day5_results_SC_SI_shrink, ylim=c(-2,2))
