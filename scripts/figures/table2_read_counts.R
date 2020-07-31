# Table of read counts.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

canola_metadata <- read.table("canola_rnaseq_metadata_edit.txt",
                              header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char="")
rownames(canola_metadata) <- canola_metadata$HS.Seq.name

fastq_lengths <- read.table("tables/fastq_line_num.txt", header=TRUE, stringsAsFactors = FALSE)
fastq_lengths$read_count <- fastq_lengths$line_numbers / 4
fastq_lengths$sample <- gsub("O.-._..*_D._", "", fastq_lengths$fastq)
fastq_lengths$sample <- gsub("_L00._R1_001.fastq.gz", "", fastq_lengths$sample)
rownames(fastq_lengths) <- fastq_lengths$sample

count_table <- fastq_lengths
count_table$aligned_reads <- NA
count_table$mmquant <- NA
count_table$mmquant_AT <- NA

flagstat_files <- list.files("bam_files/flagstat_out", full.names = TRUE)
for(flagstat in flagstat_files) {
  sample_name <- sub("bam_files/flagstat_out/", "", flagstat)
  sample_name <- sub("_flagstat.txt", "", sample_name)
  flagstat_infile <- read.delim(flagstat, header=FALSE, stringsAsFactors = FALSE, nrows=2, sep=" ")
  count_table[sample_name, "aligned_reads"] <- flagstat_infile[1, 1] - flagstat_infile[2, 1]
}


mmquant_files <- list.files("mmquant_out/mmquant_out", full.names = TRUE)
for(mmquant in mmquant_files) {
  sample_name <- sub("mmquant_out/mmquant_out/", "", mmquant)
  sample_name <- sub("_mmquant_out.txt", "", sample_name)
  mmquant_intable <- read.table(mmquant, header=TRUE, stringsAsFactors = FALSE)
  count_table[sample_name, "mmquant"] <- sum(mmquant_intable[, 2])
}

mmquant_AT_files <- list.files("mmquant_out_AT", full.names = TRUE)
for(mmquant_AT in mmquant_AT_files) {
  sample_name <- sub("mmquant_out_AT/", "", mmquant_AT)
  sample_name <- sub("_mmquant_out_AT.txt", "", sample_name)
  mmquant_AT_intable <- read.table(mmquant_AT, header=TRUE, stringsAsFactors = FALSE)
  count_table[sample_name, "mmquant_AT"] <- sum(mmquant_AT_intable[, 2])
}

# Remove unneeded columns
count_table <- count_table[, -which(colnames(count_table) %in% c("fastq", "sample", "line_numbers"))]
count_table$Sample <- canola_metadata[rownames(count_table), "Sample"]
count_table$Day <- gsub("., .., D", "D", canola_metadata$Tube.Name)

write.table(x = count_table, file = "tables/read_counts.tsv", col.names=NA, row.names=TRUE, sep="\t", quote=FALSE)

