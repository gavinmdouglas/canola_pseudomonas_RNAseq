### Commands to determine how significant genes intersect between samples and also to generate Venn diagrams.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")
source("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/scripts/canola_pseudomonas_R_code.R")

library("VennDiagram")

# Function to take in a character vector and will create output file with 1 character per line. 
# Output file name is based on input strings given in separate vector.
write_out_vec <- function(vec2write, str_for_file) {
  
  file_out <- paste(str_for_file, collapse="_")
  print(paste("Writing:", file_out))
  write.table(vec2write, file = file_out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

# Function that will draw three way Venn diagram and output gene ids overlapping in each set to textfiles.
deseq2_ThreeWayVenn_and_set <- function(results1, results2, results3, name1, name2, name3, plot_outdir,
                                        padj_cut=0.1, l2fc_cut=0, compare_type="de", 
                                        prefix = "genes", suffix="txt",
                                        venn_col=c("#1f78b4", "#33a02c", "#e31a1c")) {
  
  if(compare_type == "de") {
    # Just looking for differentially expressed. Note takes absolute of log2foldchanges.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & abs(results1$log2FoldChange) > l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & abs(results2$log2FoldChange) > l2fc_cut)]
    set3 <- rownames(results3)[which(results3$padj < padj_cut & abs(results3$log2FoldChange) > l2fc_cut)]
  } else if(compare_type == "up") {
    # Get up-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange > l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange > l2fc_cut)]
    set3 <- rownames(results3)[which(results3$padj < padj_cut & results3$log2FoldChange > l2fc_cut)]
    
  } else if(compare_type == "down") {
    # Get down-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange < l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange < l2fc_cut)]
    set3 <- rownames(results3)[which(results3$padj < padj_cut & results3$log2FoldChange < l2fc_cut)]
  } else {
    stop("Option compare_type needs to be one of \"de\", \"up\", or \"down\".")
  }
  
  # Get ids overlapping between each set.
  set_n12=set1[which(set1 %in% set2)]
  set_n23=set2[which(set2 %in% set3)]
  set_n13=set1[which(set1 %in% set3)]
  set_n123=set_n12[which(set_n12 %in% set3)]
  
  venn.plot <- draw.triple.venn(area1=length(set1),
                                area2=length(set2),
                                area3=length(set3),
                                n12=length(set_n12),
                                n23=length(set_n23),
                                n13=length(set_n13),
                                n123=length(set_n123),
                                category=c(name1, name2, name3),
                                fill=venn_col)
  
  out_info <- paste(compare_type, "padj", as.character(padj_cut), "l2fc", as.character(l2fc_cut), sep="_")
  out_suffix <- paste(out_info, suffix, sep=".")
  
  out_pdf <- paste(plot_outdir, "/", basename(prefix), "_", out_info, ".pdf", sep="")
  pdf(file=out_pdf)
  grid.draw(venn.plot);
  dev.off()
  grid.newpage();
  
  # Write out sets to files.
  write_out_vec(vec2write = set1, str_for_file = c(prefix, name1, out_suffix))
  write_out_vec(vec2write = set2, str_for_file = c(prefix, name2, out_suffix))
  write_out_vec(vec2write = set3, str_for_file = c(prefix, name3, out_suffix))
  write_out_vec(vec2write = set_n12, str_for_file = c(prefix, name1, name2, out_suffix))
  write_out_vec(vec2write = set_n23, str_for_file = c(prefix, name2, name3, out_suffix))
  write_out_vec(vec2write = set_n13, str_for_file = c(prefix, name1, name3, out_suffix))
  write_out_vec(vec2write = set_n123, str_for_file = c(prefix, name1, name2, name3, out_suffix))  
}

# Function that will draw two way Venn diagram and output gene ids overlapping in each set to textfiles.
deseq2_TwoWayVenn_and_set <- function(results1, results2, name1, name2, plot_outdir,
                                      padj_cut=0.1, l2fc_cut=0, compare_type="de", 
                                      prefix = "genes", suffix="txt",
                                      venn_col=c("#1f78b4", "#33a02c")) {
  
  if(compare_type == "de") {
    # Just looking for differentially expressed. Note takes absolute of log2foldchanges.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & abs(results1$log2FoldChange) > l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & abs(results2$log2FoldChange) > l2fc_cut)]
  } else if(compare_type == "up") {
    # Get up-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange > l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange > l2fc_cut)]
    
  } else if(compare_type == "down") {
    # Get down-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange < l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange < l2fc_cut)]
  } else {
    stop("Option compare_type needs to be one of \"de\", \"up\", or \"down\".")
  }
  
  # Get ids overlapping between each set.
  set_n12=set1[which(set1 %in% set2)]
  
  if(length(set1) < length(set2)) {
    rotation_setting <- 180
    names_set <- c(name2, name1)
  } else {
    rotation_setting <- 0
    names_set <- c(name1, name2)
  }
  
  venn.plot <- draw.pairwise.venn(
    area1=length(set1),
    area2=length(set2),
    cross.area=length(set_n12),
    category=names_set,
    cex = c(2, 2, 2),
    cat.cex = c(2, 2),
    fill=venn_col,
    rotation.degree=rotation_setting)
  
  out_info <- paste(compare_type, "padj", as.character(padj_cut), "l2fc", as.character(l2fc_cut), sep="_")
  out_suffix <- paste(out_info, suffix, sep=".")
  
  out_pdf <- paste(plot_outdir, "/", basename(prefix), "_", out_info, ".pdf", sep="")
  pdf(file=out_pdf)
  grid.draw(venn.plot);
  dev.off()
  grid.newpage();
  
  # Write out sets to files.
  write_out_vec(vec2write = set1, str_for_file = c(prefix, name1, out_suffix))
  write_out_vec(vec2write = set2, str_for_file = c(prefix, name2, out_suffix))
  write_out_vec(vec2write = set_n12, str_for_file = c(prefix, name1, name2, out_suffix))
}


dir.create("Bnapus_sig_gene_sets", showWarnings = FALSE)
dir.create("plots/Bnapus_venn", showWarnings = FALSE)


DE_categories <- c("shoot_up", "shoot_down", "shoot_de",
                   "root_up", "root_down", "root_de",
                   "day1_up", "day1_down", "day1_de",
                   "day3_up", "day3_down", "day3_de",
                   "day5_up", "day5_down", "day5_de")

for(DE_category in DE_categories) {
  dir.create(paste("Bnapus_sig_gene_sets", DE_category, sep="/"), showWarnings = FALSE)
}

# Load in DESeq2 output files.
day1_results_SC_SI_shrink <- read.table(file = "Bnapus_deseq2_outfiles/day1_results_SC_SI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)
day1_results_RC_RI_shrink <- read.table(file = "Bnapus_deseq2_outfiles/day1_results_RC_RI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)

day3_results_SC_SI_shrink <- read.table(file = "Bnapus_deseq2_outfiles/day3_results_SC_SI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)
day3_results_RC_RI_shrink <- read.table(file = "Bnapus_deseq2_outfiles/day3_results_RC_RI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)

day5_results_SC_SI_shrink <- read.table(file = "Bnapus_deseq2_outfiles/day5_results_SC_SI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)
day5_results_RC_RI_shrink <- read.table(file = "Bnapus_deseq2_outfiles/day5_results_RC_RI_shrink.txt",
                                                header=TRUE, sep="\t", row.names=1)

# Venn diagrams and gene sets.
comparison_types <- c("de", "up", "down")
l2fc_cutoffs <- 2 # c(0, 2)

for(compare in comparison_types) {
  for (l2fc in l2fc_cutoffs) {
    
    if(compare == "down") {
      l2fc = -1*l2fc
    }
    
    shoot_prefix = paste("Bnapus_sig_gene_sets/shoot_", compare, "/Bnapus_shoot_genes", sep = "")
    root_prefix = paste("Bnapus_sig_gene_sets/root_", compare, "/Bnapus_root_genes", sep = "")
    
    deseq2_ThreeWayVenn_and_set(results1 = day1_results_SC_SI_shrink,
                                results2 = day3_results_SC_SI_shrink,
                              results3 = day5_results_SC_SI_shrink,
                              name1 = "Day1",
                              name2 = "Day3",
                              name3 = "Day5",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = shoot_prefix,
                              plot_outdir="plots/Bnapus_venn")
    
    deseq2_ThreeWayVenn_and_set(results1 = day1_results_RC_RI_shrink,
                                results2 = day3_results_RC_RI_shrink,
                              results3 = day5_results_RC_RI_shrink,
                              name1 = "Day1",
                              name2 = "Day3",
                              name3 = "Day5",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = root_prefix,
                              plot_outdir="plots/Bnapus_venn")
  }
}


# Also do comparisons in between tissues.
for(compare in comparison_types) {
  for (l2fc in l2fc_cutoffs) {
    
    if(compare == "down") {
      l2fc = -1*l2fc
    }
    
    day1_prefix = paste("Bnapus_sig_gene_sets/day1_", compare, "/Bnapus_day1_genes", sep = "")
    day3_prefix = paste("Bnapus_sig_gene_sets/day3_", compare, "/Bnapus_day3_genes", sep = "")
    day5_prefix = paste("Bnapus_sig_gene_sets/day5_", compare, "/Bnapus_day5_genes", sep = "")
    
    deseq2_TwoWayVenn_and_set(results1 = day1_results_SC_SI_shrink,
                              results2 = day1_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day1_prefix,
                              plot_outdir="plots/Bnapus_venn")
    
    deseq2_TwoWayVenn_and_set(results1 = day3_results_SC_SI_shrink,
                              results2 = day3_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day3_prefix,
                              plot_outdir="plots/Bnapus_venn")
    
    deseq2_TwoWayVenn_and_set(results1 = day5_results_SC_SI_shrink,
                              results2 = day5_results_RC_RI_shrink,
                              name1 = "Shoot",
                              name2 = "Root",
                              compare_type = compare,
                              padj_cut = 0.1,
                              l2fc_cut = l2fc,
                              prefix = day5_prefix,
                              plot_outdir="plots/Bnapus_venn")
  }
}


# Count number of significant genes overall.
sig_genes <- c()

de_gene_files <- c(list.files("Bnapus_sig_gene_sets/root_de/", full.names = TRUE),
                   list.files("Bnapus_sig_gene_sets/shoot_de/", full.names = TRUE))

for(f in de_gene_files) {
  sig_genes <- c(sig_genes, read.table(f, stringsAsFactors = FALSE)$V1)
}

sig_genes <- sig_genes[-which(duplicated(sig_genes))]

write.table(x = sig_genes, file="Bnapus_sig_gene_sets/root_shoot_de.txt", col.names = FALSE, row.names=FALSE, quote=FALSE)


