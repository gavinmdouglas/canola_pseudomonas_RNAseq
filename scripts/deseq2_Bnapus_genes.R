# Commands used to run DESeq2 on RNA-seq data mapped to B. napus genes.

rm(list=ls(all.names=TRUE))

library("BiocParallel")
library("DESeq2")
library("knitr")

setwd('/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/')

# First, a sample table was prepared, which included sample information and what sequencing lane they were in.
Bn_files <- grep(".txt$", list.files("mmquant_out/mmquant_out"), value=TRUE)

samples <- gsub("_mmquant_out.txt$", "", Bn_files)

sample_table <- data.frame(sample=samples, file=Bn_files, condition=NA, stringsAsFactors=FALSE)

# Read in RNAseq mapping file as well and add condition to DF.
# My edit of this file was to convert "ul" (w special char) -> "microlitre" in header.
rnaseq_map <- read.table("canola_rnaseq_metadata_edit.txt",
                         sep="\t",
                         header=T,
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         quote="")

rownames(rnaseq_map) <- rnaseq_map$HS.Seq.name

# Get sequencing lane number for each sample
lane_map <- read.table("mapfile_lane.txt", header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE)
lane_map$lane_char <- NA
lane_map[which(lane_map$sequencing_lane == 4), "lane_char"] <- "four"
lane_map[which(lane_map$sequencing_lane == 5), "lane_char"] <- "five"

# Get corresponding condition for each sample from rnaseq_map.
sample_table$condition <- rnaseq_map[sample_table$sample, "Sample._Name"]

# Remove quotations from conditions, convert commands to hyphens,
# and remove replicate number from after C and I.
sample_table$condition <- gsub("\"", "", sample_table$condition)
sample_table$condition <- gsub(", ", "_", sample_table$condition)
sample_table$condition <- gsub("_C\\d_", "_C_", sample_table$condition)
sample_table$condition <- gsub("_I\\d_", "_I_", sample_table$condition)

# Convert all columns to be factors.
sample_table$sample <- as.factor(sample_table$sample)
sample_table$condition <- as.factor(sample_table$condition)
sample_table$file <- as.factor(sample_table$file)
sample_table$lane_char <- as.factor(lane_map[sample_table$sample, "lane_char"])


# The mmquant output was then read into a single table.
all_mmquant_categories <- c()
for(row_i in 1:nrow(sample_table)) {
  all_mmquant_categories <- c(all_mmquant_categories,
                              rownames(read.table(paste("mmquant_out/mmquant_out/", sample_table$file[row_i], sep=""), sep="\t", row.names=1, header=TRUE)))
}

all_mmquant_categories <- all_mmquant_categories[-which(duplicated(all_mmquant_categories))]

# Removed all multi-gene categories:
all_mmquant_categories <- all_mmquant_categories[-grep("--", all_mmquant_categories)]

all_mmquant_categories <- sort(all_mmquant_categories)

# Initialize empty dataframe for mmquant counts.
Bnapus_mmquant <- data.frame(matrix(0, nrow=length(all_mmquant_categories), ncol=nrow(sample_table)))
rownames(Bnapus_mmquant) <- all_mmquant_categories
colnames(Bnapus_mmquant) <- sample_table$sample

# Read in all counts.
for(row_i in 1:nrow(sample_table)) {
  sample_id <- as.character(sample_table$sample[row_i])
  mmquant_in <- read.table(paste("mmquant_out/mmquant_out/", sample_table$file[row_i], sep=""), sep="\t", row.names=1, header=TRUE)
  
  overlapping_genes <- rownames(Bnapus_mmquant)[which(rownames(Bnapus_mmquant) %in% rownames(mmquant_in))]
  Bnapus_mmquant[overlapping_genes, sample_id] <- mmquant_in[overlapping_genes, sample_id]
}


# Then ran DESeq2 based on this input table with the design ~ lane # + condition. This was done with 40 cores in parallel. The resulting objects were saved as RDS files.
sample_table_dataset <- DESeqDataSetFromMatrix(countData=Bnapus_mmquant,
                              colData = sample_table,
                              design =~ lane_char + condition)

sample_dataset <- DESeq(sample_table_dataset, parallel=TRUE, BPPARAM=MulticoreParam(40))

# Also saved main DESeq2 objects as RDS files.
saveRDS(object = sample_table_dataset, file = "Bnapus_deseq2_outfiles/sample_table.rds")
saveRDS(object = sample_dataset, file = "Bnapus_deseq2_outfiles/deseq2.rds")


#The specific comparisons per day were then conducted as well based on this output.
# Log-fold change shrinking was performed with the "ashr" method. These results for pairwise comparisons were written to tables.

sample_day1_results_SC_SI <- results(sample_dataset, contrast=c("condition", "S_I_D1", "S_C_D1"))
sample_day1_results_RC_RI <- results(sample_dataset, contrast=c("condition", "R_I_D1", "R_C_D1"))

sample_day3_results_SC_SI <- results(sample_dataset, contrast=c("condition", "S_I_D3", "S_C_D3"))
sample_day3_results_RC_RI <- results(sample_dataset, contrast=c("condition", "R_I_D3", "R_C_D3"))

sample_day5_results_SC_SI <- results(sample_dataset, contrast=c("condition", "S_I_D5", "S_C_D5"))
sample_day5_results_RC_RI <- results(sample_dataset, contrast=c("condition", "R_I_D5", "R_C_D5"))



# Shrink log-fold change values.
day1_results_SC_SI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "S_I_D1", "S_C_D1"), 
                                               res=sample_day1_results_SC_SI, type="ashr")

day1_results_RC_RI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "R_I_D1", "R_C_D1"), 
                                               res=sample_day1_results_RC_RI, type="ashr")

day3_results_SC_SI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "S_I_D3", "S_C_D3"), 
                                               res=sample_day3_results_SC_SI, type="ashr")

day3_results_RC_RI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "R_I_D3", "R_C_D3"), 
                                               res=sample_day3_results_RC_RI, type="ashr")

day5_results_SC_SI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "S_I_D5", "S_C_D5"), 
                                               res=sample_day5_results_SC_SI, type="ashr")

day5_results_RC_RI_shrink <- lfcShrink(sample_dataset, 
                                               contrast=c("condition", "R_I_D5", "R_C_D5"), 
                                               res=sample_day5_results_RC_RI, type="ashr")


### Write out intermediate files after shrinking l2fc values.
write.table(day1_results_SC_SI_shrink, file="Bnapus_deseq2_outfiles/day1_results_SC_SI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day1_results_RC_RI_shrink, file="Bnapus_deseq2_outfiles/day1_results_RC_RI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(day3_results_SC_SI_shrink, file="Bnapus_deseq2_outfiles/day3_results_SC_SI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day3_results_RC_RI_shrink, file="Bnapus_deseq2_outfiles/day3_results_RC_RI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(day5_results_SC_SI_shrink, file="Bnapus_deseq2_outfiles/day5_results_SC_SI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(day5_results_RC_RI_shrink, file="Bnapus_deseq2_outfiles/day5_results_RC_RI_shrink.txt", 
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)



# Create table with the log-fold change values for all pairwise comparisons and add in A. thaliana annotation information when available.
# This table was written to an output file.

combined_deseq2_out <- cbind(        day1_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                                     day3_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                                     day5_results_SC_SI_shrink[,"log2FoldChange", drop=FALSE],
                                     day1_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE],
                                     day3_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE],
                                     day5_results_RC_RI_shrink[,"log2FoldChange", drop=FALSE])

colnames(combined_deseq2_out) <- c("s1", "s3", "s5", "r1", "r3", "r5")

Bnapus_annot_map <- read.table("tables/Bnapus_merged_func_annot.txt",
                               sep="\t",
                               header=T,
                               stringsAsFactors = FALSE,
                               comment.char = "",
                               quote="")

rownames(Bnapus_annot_map) <- make.names(Bnapus_annot_map$Bnapus_gene_symbol, unique=TRUE)

combined_deseq2_out$Bnapus_genes <- rownames(combined_deseq2_out)
combined_deseq2_out$Athaliana_top_hit <- Bnapus_annot_map[rownames(combined_deseq2_out), "Athaliana_gene"]
combined_deseq2_out$Athaliana_description <- Bnapus_annot_map[rownames(combined_deseq2_out), "Athaliana_description"]

combined_deseq2_out <- combined_deseq2_out[, c("Bnapus_genes", "Athaliana_top_hit", "Athaliana_description", "s1", "s3", "s5", "r1", "r3", "r5")]

write.table(combined_deseq2_out, file="Bnapus_deseq2_outfiles/deseq2_log2fold.txt", 
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

