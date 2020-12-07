# Figure of log2 fold difference for key immune regulator genes

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq")

library(cowplot)
library(ggbeeswarm)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

Bnapus_annot <- read.table("tables/Bnapus_merged_func_annot.txt",
                           header=TRUE, sep="\t", row.names=1, quote="", comment.char="", stringsAsFactors = FALSE)


# Get B napus genes that match At homolog. If more than two matches then take the top two with the highest ID matches.
# Bnapus_BLAST <- read.table("tables/blastn_out_napus_vs_thaliana_evalue0.0001_clean_tophits.txt",
#                            header=TRUE, sep="\t", row.names=1, quote="", comment.char="", stringsAsFactors = FALSE)
# Bnapus_BLAST$Bn_gene <- rownames(Bnapus_BLAST)
# rownames(Bnapus_BLAST) <- Bnapus_annot[Bnapus_BLAST$Bn_gene, "Bnapus_gene_symbol"]
# NOTE: In the end decided to just show all homologues anyway.
# at_most_two_top_homologs <- function(gene_matches) {
#   if(length(gene_matches) > 2) {
#     BLAST_subset <- Bnapus_BLAST[gene_matches, ]
#     gene_matches <- rownames(BLAST_subset[order(BLAST_subset$pident, decreasing = TRUE),])[c(1, 2)]
#   }
#   
#   return(gene_matches)
# }

# Get Bnapus gene sets for each At gene of interest (innate immunity TFs)

gene_to_ids <- list()

# ERF094
gene_to_ids[["ERF094"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "ERF094"), "Bnapus_gene_symbol"]

# ERF1
gene_to_ids[["ERF1-1"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "ERF1-1"), "Bnapus_gene_symbol"]

# ERF6
gene_to_ids[["ERF6"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "ERF6"), "Bnapus_gene_symbol"]

# ERF104
gene_to_ids[["ERF104"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "ERF104"), "Bnapus_gene_symbol"]

# MYC2
gene_to_ids[["MYC2"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "MYC2"), "Bnapus_gene_symbol"]

# MYC3
gene_to_ids[["MYC3"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "MYC3"), "Bnapus_gene_symbol"]

# MYC4
gene_to_ids[["MYC4"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "MYC4"), "Bnapus_gene_symbol"]


# TGA1 and TGA4 (clade I)
gene_to_ids[["TGA1"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "TGA1"), "Bnapus_gene_symbol"]
gene_to_ids[["TGA4"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "TGA4"), "Bnapus_gene_symbol"]

# TGA2, TGA5, and TGA6 (clade II)
gene_to_ids[["TGA2"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "TGA2"), "Bnapus_gene_symbol"]
gene_to_ids[["TGA5"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "TGA5"), "Bnapus_gene_symbol"]
gene_to_ids[["TGA6"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "TGA6"), "Bnapus_gene_symbol"]

# TGA3 and TGA7 (clade III)
gene_to_ids[["TGA3"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "TGA3"), "Bnapus_gene_symbol"]
gene_to_ids[["TGA7"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "TGA7"), "Bnapus_gene_symbol"]

# MYB family
gene_to_ids[["MYB30"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "MYB30"), "Bnapus_gene_symbol"]
gene_to_ids[["MYB44"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "MYB44"), "Bnapus_gene_symbol"]
gene_to_ids[["MYB108"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "MYB108"), "Bnapus_gene_symbol"]

# NAC family
gene_to_ids[["NAC019"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "NAC019"), "Bnapus_gene_symbol"]
gene_to_ids[["NAC055"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "NAC055"), "Bnapus_gene_symbol"]
# Couldn't find any matches to NAC072

# WRKY family
gene_to_ids[["WRKY22"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "WRKY22"), "Bnapus_gene_symbol"]
gene_to_ids[["WRKY29"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "WRKY29"), "Bnapus_gene_symbol"]
gene_to_ids[["WRKY30"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "WRKY30"), "Bnapus_gene_symbol"]
gene_to_ids[["WRKY33"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "WRKY33"), "Bnapus_gene_symbol"]
# Note that no hits for WRKY52 identified.
#gene_to_ids[["WRKY52"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "WRKY52"), "Bnapus_gene_symbol"]
gene_to_ids[["WRKY53"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "WRKY53"), "Bnapus_gene_symbol"]
gene_to_ids[["WRKY70"]] <- Bnapus_annot[which(Bnapus_annot$Athaliana_gene_symbol == "WRKY70"), "Bnapus_gene_symbol"]


all_Bnapus_genes <- c()
for(gene_id in names(gene_to_ids)) {
  
  # Could use this command to get two best hits, but decided to just show all possible homologues.
  # if(length(gene_to_ids[[gene_id]]) > 2) {
  #   gene_to_ids[[gene_id]] <- at_most_two_top_homologs(gene_to_ids[[gene_id]])
  # }
  
  all_Bnapus_genes <- c(all_Bnapus_genes, gene_to_ids[[gene_id]])
}


# Double-check a few At gene id matches and that they are the top BLAST hits if > 2.
# ERF104 hits: c("BnaA03g40380D", "BnaA08g13940D", "BnaC07g31350D")
# Bnapus_annot[which(Bnapus_annot$Bnapus_gene_symbol %in% c("BnaA03g40380D", "BnaA08g13940D", "BnaC07g31350D")), ]
# Bnapus_BLAST[c("BnaA03g40380D", "BnaA08g13940D", "BnaC07g31350D"), ]
# gene_to_ids[["ERF104"]]
# 
# # TGA3 hits: c("BnaA08g20970D", "BnaA09g30910D", "BnaC05g17700D", "BnaC08g20170D")
# Bnapus_annot[which(Bnapus_annot$Bnapus_gene_symbol %in% c("BnaA08g20970D", "BnaA09g30910D", "BnaC05g17700D", "BnaC08g20170D")), ]
# Bnapus_BLAST[c("BnaA08g20970D", "BnaA09g30910D", "BnaC05g17700D", "BnaC08g20170D"), ]
# gene_to_ids[["TGA3"]]
# 
# # WRKY29 hits: c("BnaA01g13230D", "BnaA03g46020D", "BnaA08g10660D", "BnaC01g15320D", "BnaC03g65200D", "BnaC07g25890D", "BnaC07g38350D")
# Bnapus_annot[which(Bnapus_annot$Bnapus_gene_symbol %in% c("BnaA01g13230D", "BnaA03g46020D", "BnaA08g10660D", "BnaC01g15320D", "BnaC03g65200D", "BnaC07g25890D", "BnaC07g38350D")), ]
# Bnapus_BLAST[c("BnaA01g13230D", "BnaA03g46020D", "BnaA08g10660D", "BnaC01g15320D", "BnaC03g65200D", "BnaC07g25890D", "BnaC07g38350D"), ]
# gene_to_ids[["WRKY29"]]


Bnapus_deseq2_out <- read.table("Bnapus_deseq2_outfiles/deseq2_log2fold.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", row.names=1)

rna_levels <- list()

for(tissue in c("root", "shoot")) {
  rna_levels[[tissue]] <- list()
  
  for(day in c("day1", "day3", "day5")) {
    
    rna_levels[[tissue]][[day]] <- data.frame(matrix(NA, nrow=length(all_Bnapus_genes), ncol=6))
    colnames(rna_levels[[tissue]][[day]]) <- c("Bn_gene", "At_gene", "Day", "log2fold", "log2fold_se", "padj")
    rownames(rna_levels[[tissue]][[day]]) <- all_Bnapus_genes
    rna_levels[[tissue]][[day]]$Bn_gene <- all_Bnapus_genes
    rna_levels[[tissue]][[day]]$Day <- day
    
    if(tissue == "root") {
      lfcshrink_file <- paste("Bnapus_deseq2_outfiles/", day, "_results_RC_RI_shrink.txt", sep="")
    } else if(tissue == "shoot") {
      lfcshrink_file <- paste("Bnapus_deseq2_outfiles/", day, "_results_SC_SI_shrink.txt", sep="")
    }

    tmp_lfcshrink <- read.table(lfcshrink_file, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
    
    for(At_gene in names(gene_to_ids)) {
      rna_levels[[tissue]][[day]][gene_to_ids[[At_gene]], "At_gene"] <- At_gene
      rna_levels[[tissue]][[day]][gene_to_ids[[At_gene]], "log2fold"] <- tmp_lfcshrink[gene_to_ids[[At_gene]], "log2FoldChange"]
      rna_levels[[tissue]][[day]][gene_to_ids[[At_gene]], "log2fold_se"] <- tmp_lfcshrink[gene_to_ids[[At_gene]], "lfcSE"]
      rna_levels[[tissue]][[day]][gene_to_ids[[At_gene]], "padj"] <- tmp_lfcshrink[gene_to_ids[[At_gene]], "padj"]
    }
  }
}

root_rna <- do.call("rbind", rna_levels$root)
shoot_rna <- do.call("rbind", rna_levels$shoot)

root_rna[which(root_rna$Day == "day1"), "Day"] <- "1"
root_rna[which(root_rna$Day == "day3"), "Day"] <- "3"
root_rna[which(root_rna$Day == "day5"), "Day"] <- "5"

root_rna$padj_category <- ">= 0.1"
root_rna[which(root_rna$padj < 0.1), "padj_category"] <- "< 0.1"
root_rna[which(root_rna$padj < 0.01), "padj_category"] <- "< 0.01"

shoot_rna[which(shoot_rna$Day == "day1"), "Day"] <- "1"
shoot_rna[which(shoot_rna$Day == "day3"), "Day"] <- "3"
shoot_rna[which(shoot_rna$Day == "day5"), "Day"] <- "5"

shoot_rna$padj_category <- ">= 0.1"
shoot_rna[which(shoot_rna$padj < 0.1), "padj_category"] <- "< 0.1"
shoot_rna[which(shoot_rna$padj < 0.01), "padj_category"] <- "< 0.01"


high_effect_genes <- c(root_rna[which(abs(root_rna$log2fold) > 1.5), "At_gene"],
                       shoot_rna[which(abs(shoot_rna$log2fold) > 1.5), "At_gene"])
high_effect_genes <- high_effect_genes[-which(duplicated(high_effect_genes))]

root_rna_high_effect <- root_rna[which(root_rna$At_gene %in% high_effect_genes), ]
shoot_rna_high_effect <- shoot_rna[which(shoot_rna$At_gene %in% high_effect_genes), ]

root_rna_plot <- ggplot(root_rna, aes(x=Day, y=log2fold, fill=padj_category)) + 
                         geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.5) +
                           geom_quasirandom(size=1.5, pch=21) +
                         facet_wrap(~ At_gene, nrow = 4) +
                         ggtitle("Root") +
                         xlab("Day") +
                         ylab(expression('log'[2]*' fold change')) +
                         scale_fill_manual(name="Adjusted\nP-value", values=c("black", "grey75", "white")) +
                         theme_bw() +
                         theme(panel.grid.minor = element_blank(),
                               plot.title = element_text(hjust = 0.5)) +
                    ylim(c(-4, 8))


root_rna_high_effect_plot <- ggplot(root_rna_high_effect, aes(x=Day, y=log2fold, fill=padj_category)) + 
                                    geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.5) +
                                    geom_quasirandom(size=1.5, pch=21) +
                                    facet_wrap(~ At_gene, nrow = 2) +
                                    ggtitle("Root") +
                                    xlab("Day") +
                                    ylab(expression('log'[2]*' fold change')) +
                                    scale_fill_manual(name="Adjusted\nP-value", values=c("black", "grey75", "white")) +
                                    theme_bw() +
                                    theme(panel.grid.minor = element_blank()) +
                                    ylim(c(-4, 8))

shoot_rna_plot <- ggplot(shoot_rna, aes(x=Day, y=log2fold, fill=padj_category)) + 
                          geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.5) +
                          geom_quasirandom(size=1.5, pch=21) +
                          facet_wrap(~ At_gene, nrow = 4) +
                          ggtitle("Shoot") +
                          xlab("Day") +
                          ylab(expression('log'[2]*' fold change')) +
                          scale_fill_manual(name="Adjusted\nP-value", values=c("black", "grey75", "white")) +
                          theme_bw() +
                          theme(panel.grid.minor = element_blank(),
                                plot.title = element_text(hjust = 0.5)) +
                          ylim(c(-4, 8))

shoot_rna_high_effect_plot <- ggplot(shoot_rna_high_effect, aes(x=Day, y=log2fold, fill=padj_category)) + 
                                    geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.5) +
                                    geom_quasirandom(size=1.5, pch=21) +
                                    facet_wrap(~ At_gene, nrow = 2) +
                                    ggtitle("Shoot") +
                                    xlab("Day") +
                                    ylab(expression('log'[2]*' fold change')) +
                                    scale_fill_manual(name="Adjusted\nP-value", values=c("black", "grey75",  "white")) +
                                    theme_bw() +
                                    theme(panel.grid.minor = element_blank()) +
                                    ylim(c(-4, 8))

combined_high_effect_plot <- plot_grid(root_rna_high_effect_plot, shoot_rna_high_effect_plot, labels=c('A', 'B'), nrow=2)

ggsave(filename = "plots/main/Figure8_master_TFs.pdf", plot = combined_high_effect_plot,
       width = 7.20472, height=5)

tiff(file = "plots/main/Figure8_master_TFs.tiff", width=7.20472, height=5,
     compression = "lzw", res=300, units="in")
combined_high_effect_plot
dev.off()

ggsave(filename = "plots/exploratory/all_shoot_master_TFs.pdf", plot = shoot_rna_plot,
       width = 7.20472, height=5)

ggsave(filename = "plots/exploratory/all_root_master_TFs.pdf", plot = root_rna_plot,
       width = 7.20472, height=5)

