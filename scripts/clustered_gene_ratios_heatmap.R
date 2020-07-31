# This script was used to read in the annotated log fold ratios, cluster the 
# genes, and to plot a heatmap and dendogram. Note that only genes with 
# A. thaliana homologs will be retained. Output files were also generated that 
# can be used with cluster 3.

#Package to explore missingness in data.
library(Amelia)

# Package for heatmap.2.
library(gplots)

# Package for splitting 1 column into multiple.
library(stringr)

# Set working directory.
setwd("Dropbox/work/Langille/pseudomonas/canola_pseudomonas/")

# Read in and combine all original differentially expressed (DE) gene sets.
root1 <- read.table("RNAseq_DE_genes/original_output/RootDay1.txt", 
                    sep="\t", stringsAsFactors = FALSE, comment.char = "", 
                    quote="", header=TRUE)

root3 <- read.table("RNAseq_DE_genes/original_output/RootDay3.txt", 
                    sep="\t", stringsAsFactors = FALSE, comment.char = "", 
                    quote="", header=TRUE)

root5 <- read.table("RNAseq_DE_genes/original_output/RootDay5.txt", 
                    sep="\t", stringsAsFactors = FALSE, comment.char = "", 
                    quote="", header=TRUE)


shoot1 <- read.table("RNAseq_DE_genes/original_output/ShootDay1.txt", 
                    sep="\t", stringsAsFactors = FALSE, comment.char = "", 
                    quote="", header=TRUE)

shoot3 <- read.table("RNAseq_DE_genes/original_output/ShootDay3.txt", 
                    sep="\t", stringsAsFactors = FALSE, comment.char = "", 
                    quote="", header=TRUE)

shoot5 <- read.table("RNAseq_DE_genes/original_output/ShootDay5.txt", 
                    sep="\t", stringsAsFactors = FALSE, comment.char = "", 
                    quote="", header=TRUE)

# Remove p-value from each df.
root1 <- root1[,-3]
root3 <- root3[,-3]
root5 <- root5[,-3]
shoot1 <- shoot1[,-3]
shoot3 <- shoot3[,-3]
shoot5 <- shoot5[,-3]

# Change gene and log2foldchange column names to be consistent.
colnames(root1) <- c("gene", "root1_log_fold")
colnames(root3) <- c("gene", "root3_log_fold")
colnames(root5) <- c("gene", "root5_log_fold")
colnames(shoot1) <- c("gene", "shoot1_log_fold")
colnames(shoot3) <- c("gene", "shoot3_log_fold")
colnames(shoot5) <- c("gene", "shoot5_log_fold")

# Put all of these dataframes (dfs) into a list so that Reduce function can be 
# used below.
DE_genes_list <- list("root1"=root1, 
                      "root3"=root3, 
                      "root5"=root5,
                      "shoot1"=shoot1, 
                      "shoot3"=shoot3, 
                      "shoot5"=shoot5)

# Merge all tables into 1 df.
merged_DE_genes <- Reduce(function(x, y) merge(x, y, all=TRUE, 
                          by="gene"), DE_genes_list, accumulate=FALSE)

# Read in combined annotation too.
combined_annotation <- read.table("tables/Bnapus_merged_func_annot.txt",
                                  header=T, stringsAsFactors = F, sep="\t", 
                                  comment.char="", quote="")

# Restrict to gene name, At homolog, and homolog description.
combined_annotation <- combined_annotation[,c("Bnapus_gene", "Athaliana_gene", 
                                              "Athaliana_description")]
merged_DE_genes_annot <- merge(merged_DE_genes,
                               combined_annotation,
                               by.x="gene",
                               by.y="Bnapus_gene",
                               all=TRUE)

# Exclude all rows that don't have an A. thaliana homolog.
merged_DE_genes_annot_At <- merged_DE_genes_annot[-which(is.na(merged_DE_genes_annot$Athaliana_gene)),]

# Make new column of combined At gene and description.
merged_DE_genes_annot_At$At_info <- as.factor(paste(merged_DE_genes_annot_At$Athaliana_gene, 
                                          merged_DE_genes_annot_At$Athaliana_description,
                                          sep="_"))
# Remove other metadata columns.
merged_DE_genes_annot_At <- merged_DE_genes_annot_At[, -c(1,8,9)]

# There are so many NA values in this dataset, which will cause problems for
# clustering. The NA values across the df can be viewed with this command:

# missmap(merged_DE_genes_annot_At)

# Replace all missing values with 0. This is a naive way to fix the issue.
merged_DE_genes_annot_At[is.na(merged_DE_genes_annot_At)] <- 0

# Get mean log ratio for all genes that share the same A. thaliana homolog.
merged_DE_genes_annot_At_mean <- aggregate(.~At_info, 
                                           data=merged_DE_genes_annot_At,
                                           mean)
# Make output file for cluster3.
merged_DE_genes_annot_At_mean_cluster3 <- merged_DE_genes_annot_At_mean

separated_AT_info <- data.frame(str_split_fixed(merged_DE_genes_annot_At_mean_cluster3$At_info, "_", 2))
colnames(separated_AT_info) <- c("At_gene", "NAME")

merged_DE_genes_annot_At_mean_cluster3 <- merged_DE_genes_annot_At_mean_cluster3[,-1]
merged_DE_genes_annot_At_mean_cluster3 <- cbind(separated_AT_info,
                                                merged_DE_genes_annot_At_mean_cluster3)

write.table(x = merged_DE_genes_annot_At_mean_cluster3, file = "tables/merged_DE_cluster3_input.txt", col.names=TRUE, row.names=FALSE, quote = FALSE, sep="\t")

# Set At_info to be rownames and remove column.
rownames(merged_DE_genes_annot_At_mean) <- merged_DE_genes_annot_At_mean$At_info
merged_DE_genes_annot_At_mean <- merged_DE_genes_annot_At_mean[, -1]

# Simplify column names for plotting
colnames(merged_DE_genes_annot_At_mean) <- gsub("_log_fold", "", 
                                                colnames(merged_DE_genes_annot_At_mean))

merged_DE_genes_annot_At_mean_dist <- dist(merged_DE_genes_annot_At_mean)

merged_DE_genes_annot_At_mean_dist_hclust <- hclust(merged_DE_genes_annot_At_mean_dist)

merged_DE_genes_annot_At_mean_matrix <- as.matrix(merged_DE_genes_annot_At_mean)

heatmap.2(as.matrix(merged_DE_genes_annot_At_mean),
          dendrogram="none",trace="none", scale="none")

#heatmap.2(as.matrix(merged_DE_genes_annot_At_mean), col=redgreen(75), 
#          density.info="none", trace="none", 
#          Rowv = as.dendrogram(merged_DE_genes_annot_At_mean_dist_hclust), 
#          symm=F,symkey=T,symbreaks=T, scale="none") 

heatmap.2(as.matrix(merged_DE_genes_annot_At_mean),
          margin=c(5,15),dendrogram="none",trace="none",scale="none",
          Rowv=FALSE)



          
