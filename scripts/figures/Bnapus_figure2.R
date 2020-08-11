# Figure of summary heatmap

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

library(circlize)
library(ComplexHeatmap)
library(Hmisc)

# Read in log2fold ratios.
ratios_Bn <- read.table("Bnapus_deseq2_outfiles/deseq2_log2fold.txt",
                        header=TRUE, 
                        sep="\t",
                        stringsAsFactors = FALSE,
                        quote="",
                        comment.char = "")

rownames(ratios_Bn) <- ratios_Bn$Bnapus_genes

# Get info columns and remove from original df.
ratios_Bn_info <- ratios_Bn[, c(1, 2, 3)]
ratios_Bn <- ratios_Bn[, -c(1, 2, 3)]

# Make all descriptions unique.
ratios_Bn_info$Athaliana_description_unique <- make.unique(ratios_Bn_info$Athaliana_description, sep = ".")

# Read in significant genes in at least 1 tissue on at least 1 day.
Bnapus_sig_genes <- read.table("Bnapus_sig_gene_sets/root_shoot_de.txt", header=FALSE, stringsAsFactors = FALSE)$V1

# Subset to the rows idnetified above and set any absolute values greater than 2 to have a max abs value of 2 (to make it easier to visualize).
ratios_Bn_set <- ratios_Bn[Bnapus_sig_genes, ]
ratios_Bn_info <- ratios_Bn_info[Bnapus_sig_genes, ]

ratios_Bn_set[ratios_Bn_set > 2] <- 2
ratios_Bn_set[ratios_Bn_set < -2] <- -2

# Re-order so that root samples are first.
ratios_Bn_set <- ratios_Bn_set[, c("r1", "r3", "r5", "s1", "s3", "s5")]

# Give columns clearer names
colnames(ratios_Bn_set) <- c("Day 1 Root", "Day 3 Root", "Day 5 Root", "Day 1 Shoot", "Day 3 Shoot", "Day 5 Shoot")

# Read in link from At genes to biological processes.
At_gene_single_category <- readRDS("At_GO_RDS/At_gene_single_category.rds")

# First figure out matches of all B. napus hits to each category.
# Any categories with fewer than 100 genes will be collapsed to the "Other biological process" categeory.
# Genes in these categories will be mapped to the next smallest functional category if possible or will be retained in the "Other biological process" category.

Bnapus_At_categories <- list()

for(category in names(At_gene_single_category)) {
  Bnapus_At_categories[[category]] <- character()
}


matched_Bnapus_genes <- rownames(ratios_Bn_info)[-which(is.na(ratios_Bn_info$Athaliana_top_hit))]

for(matched_gene in matched_Bnapus_genes) {
  At_gene <- ratios_Bn_info[matched_gene, "Athaliana_top_hit"]

  for(category in names(sort(sapply(At_gene_single_category, length), decreasing = TRUE))) {
     if(At_gene %in% At_gene_single_category[[category]]) {
       Bnapus_At_categories[[category]] <- c(Bnapus_At_categories[[category]], At_gene)
       break
     }
  }

}

Bnapus_At_category_lengths <- sapply(Bnapus_At_categories, length)

Bnapus_At_category_other <- names(Bnapus_At_category_lengths)[which(Bnapus_At_category_lengths < 100)]

At_genes_to_GO <- read.table("At_GO/ATH_to_GO_unique.tsv",
                             header=FALSE, sep="\t", stringsAsFactors = FALSE)

GO_slim_to_At <- readRDS("At_GO_RDS/GO_slim_to_At.rds")

group_At_genes_to_key_BP <- function(At_genes_to_GO,
                                     GO_slim_to_At,
                                     other_cat=c()) {

  category_lengths <- sort(sapply(GO_slim_to_At, length))
  
  At_genes <- At_genes_to_GO$V1[-which(duplicated(At_genes_to_GO$V1))]
  
  At_gene_single_category <- list()
  
  for(category in names(category_lengths)) {
    if(category %in% other_cat) { next }
    At_gene_single_category[[category]] <- as.character()
  }
  
  At_gene_single_category[["Unknown biological processes"]] <- as.character()
  At_gene_single_category[["Other biological processes"]] <- as.character()
  
  for(At in At_genes) {
    
    At_categorized <- FALSE
    other_category <- FALSE
    
    for(category in names(category_lengths)) {
      if(At %in% names(GO_slim_to_At[[category]])) {
        
        if(category %in% other_cat) {
          other_category <- TRUE
          next
        }
        
        At_gene_single_category[[category]] <- c(At_gene_single_category[[category]], At)
        At_categorized <- TRUE
        break
      }
    }
    
    if(! At_categorized) {
      if(other_category) {
        At_gene_single_category[["Other biological processes"]] <- c(At_gene_single_category[["Other biological processes"]], At)
      } else {
        At_gene_single_category[["Unknown biological processes"]] <- c(At_gene_single_category[["Unknown biological processes"]], At)
      }
    }
  }
  
  return(At_gene_single_category)
}


# Read in GO slim table.
GO_main_categories <- read.table("At_GO/TAIR_GO_slim_categories.txt",
                                 comment.char = "!", sep="\t", skip=8, quote="", stringsAsFactors = FALSE, header=TRUE)
# Remove quote characters.
for(coln in colnames(GO_main_categories)) {
  GO_main_categories[, coln] <- gsub("\"", "", GO_main_categories[, coln])
}
GO_main_bp <- GO_main_categories[which(GO_main_categories$ONTOLOGY.ASPECT == "biological process"), ]
rownames(GO_main_bp) <- GO_main_bp$GO_ID

At_gene_single_category_collapsed <- group_At_genes_to_key_BP(At_genes_to_GO = At_genes_to_GO,
                                                              GO_slim_to_At = GO_slim_to_At,
                                                              other_cat = Bnapus_At_category_other)

split_vec <- c()
for(Bn_gene in rownames(ratios_Bn_set)) {

  At_gene <- ratios_Bn_info[Bn_gene, "Athaliana_top_hit"]
  
  annotated <- FALSE

  for(category in names(sort(sapply(At_gene_single_category_collapsed, length), decreasing = TRUE))) {

    if(At_gene %in% At_gene_single_category_collapsed[[category]]) {
      
      if(category %in% rownames(GO_main_bp)) {
        category_descrip <- GO_main_bp[category, "SLIM_NAME"]
      } else {
        category_descrip <- category
      }
      
      # Merge unknown biological processes with "Other".
      if(category_descrip == "Unknown biological processes") {
        category_descrip <- "Other biological processes"
      }

      split_vec <- c(split_vec, capitalize(category_descrip))
      annotated <- TRUE
      break
    }
  }
  
  if(! annotated) {
    split_vec <- c(split_vec, "No BLAST match")  
  }
}

split_vec[which(split_vec == "Nucleobase-containing compound metabolic process")] <- "Nucleobase-containing comp. met. proc."

col_rnorm = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

full_heatmap <- Heatmap(matrix=as.matrix(ratios_Bn_set),
                        col=col_rnorm,
                        name="Log fold change",
                        cluster_columns = FALSE,
                        show_row_dend=TRUE,
                        show_row_names=FALSE,
                        split=split_vec,
                        column_split=c(rep("Root", 3), rep("Shoot", 3)),
                        column_names_rot = 45,
                        row_title_rot=0,
                        row_gap=unit(3, "mm"),
                        row_title_side="right",
                        column_gap=unit(4, "mm"),
                        clustering_distance_rows="euclidean")

pdf(file = "plots/main/Figure2.pdf", width=7.3, height=7.3, onefile=FALSE)
full_heatmap
dev.off()
