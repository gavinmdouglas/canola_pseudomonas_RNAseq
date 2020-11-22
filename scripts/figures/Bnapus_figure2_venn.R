# Figure of main venn diagrams

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

library(cowplot)
library(ggVennDiagram)


# Venn diagrams for each tissue of overall DE genes (lfc > 2) by day. Four panels in total for up/down in both shoots and roots

# Panel A - root up Venn diagram
root_up_sig <- list()
root_up_sig[["Day 1"]] <- read.table("Bnapus_sig_gene_sets/day1_up/Bnapus_day1_genes_Root_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
root_up_sig[["Day 3"]] <- read.table("Bnapus_sig_gene_sets/day3_up/Bnapus_day3_genes_Root_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
root_up_sig[["Day 5"]] <- read.table("Bnapus_sig_gene_sets/day5_up/Bnapus_day5_genes_Root_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

root_up_sig_venn <- ggVennDiagram(x = root_up_sig) +
                        scale_fill_gradient(low="light grey", high = "red", limits = c(0, 2000)) +
                        labs(fill="Count") +
                        ggtitle("Root up-regulated") +
  theme(plot.title = element_text(hjust = 0.5))

# Panel B - root down Venn diagram
root_down_sig <- list()
root_down_sig[["Day 1"]] <- read.table("Bnapus_sig_gene_sets/day1_down/Bnapus_day1_genes_Root_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1
root_down_sig[["Day 3"]] <- read.table("Bnapus_sig_gene_sets/day3_down/Bnapus_day3_genes_Root_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1
root_down_sig[["Day 5"]] <- read.table("Bnapus_sig_gene_sets/day5_down/Bnapus_day5_genes_Root_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

root_down_sig_venn <- ggVennDiagram(x = root_down_sig) +
                          scale_fill_gradient(low="light grey", high = "red", limits = c(0, 2000)) +
                          labs(fill="Count") +
                          ggtitle("Root down-regulated") +
  theme(plot.title = element_text(hjust = 0.5))


# Panel C - shoot up Venn diagram
shoot_up_sig <- list()
shoot_up_sig[["Day 1"]] <- read.table("Bnapus_sig_gene_sets/day1_up/Bnapus_day1_genes_Shoot_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
shoot_up_sig[["Day 3"]] <- read.table("Bnapus_sig_gene_sets/day3_up/Bnapus_day3_genes_Shoot_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1
shoot_up_sig[["Day 5"]] <- read.table("Bnapus_sig_gene_sets/day5_up/Bnapus_day5_genes_Shoot_up_padj_0.1_l2fc_2.txt", stringsAsFactors = FALSE)$V1

shoot_up_sig_venn <- ggVennDiagram(x = shoot_up_sig) +
  scale_fill_gradient(low="light grey", high = "red", limits = c(0, 2000)) +
  labs(fill="Count") +
  ggtitle("Shoot up-regulated") +
  theme(plot.title = element_text(hjust = 0.5))

# Panel D - shoot down Venn diagram
shoot_down_sig <- list()
shoot_down_sig[["Day 1"]] <- read.table("Bnapus_sig_gene_sets/day1_down/Bnapus_day1_genes_Shoot_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1
shoot_down_sig[["Day 3"]] <- read.table("Bnapus_sig_gene_sets/day3_down/Bnapus_day3_genes_Shoot_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1
shoot_down_sig[["Day 5"]] <- read.table("Bnapus_sig_gene_sets/day5_down/Bnapus_day5_genes_Shoot_down_padj_0.1_l2fc_-2.txt", stringsAsFactors = FALSE)$V1

shoot_down_sig_venn <- ggVennDiagram(x = shoot_down_sig) +
  scale_fill_gradient(low="light grey", high = "red", limits = c(0, 2000)) +
  labs(fill="Count") +
  ggtitle("Shoot down-regulated") +
  theme(plot.title = element_text(hjust = 0.5))



# Plot figure
pdf(file = "plots/main/Figure2_venn.pdf", width=7.3, height=7.3, onefile=FALSE)
plot_grid(root_up_sig_venn, root_down_sig_venn, shoot_up_sig_venn, shoot_down_sig_venn,
          labels=c('A', 'B', 'C', 'D'), nrow=2, ncol=2)
dev.off()


tiff(file = "plots/main/Figure2_venn.tiff", width=7.3, height=7.3,
     compression = "lzw", res=300, units="in")
plot_grid(root_up_sig_venn, root_down_sig_venn, shoot_up_sig_venn, shoot_down_sig_venn,
          labels=c('A', 'B', 'C', 'D'), nrow=2, ncol=2)
dev.off()
