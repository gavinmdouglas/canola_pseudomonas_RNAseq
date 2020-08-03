rm(list=ls(all.names=TRUE))

library("GOfuncR")
library("knitr")

setwd('/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/')

# Categorize all A. thaliana genes into high-level GO biological processes.
# These categorizes are used by TAIR for large-scale visualization. The links from these categories to genes will be output as a list in RDS format.
# A simplified output will also be made where each A. thaliana gene is only included in a single category.
# This will be done by considering genes as only part of the smallest category that it is part of (which greatly simplifies how the data can be displayed).
# This will also be saved as an RDS file.

# First read in all key GO categories that will be visualized (limited to biological processes).
GO_main_categories <- read.table("At_GO/TAIR_GO_slim_categories.txt",
                                 comment.char = "!", sep="\t", skip=8, quote="", stringsAsFactors = FALSE, header=TRUE)

# Remove quote characters.
for(coln in colnames(GO_main_categories)) {
  GO_main_categories[, coln] <- gsub("\"", "", GO_main_categories[, coln])
}

GO_main_bp <- GO_main_categories[which(GO_main_categories$ONTOLOGY.ASPECT == "biological process"), ]

rownames(GO_main_bp) <- GO_main_bp$GO_ID

# Acquired all GO id subcategories within each of these GO SLIM categories with the GOfuncR package.

go_slim_subcategories <- lapply(GO_main_bp$GO_ID, get_child_nodes)
names(go_slim_subcategories) <- GO_main_bp$GO_ID


# Read in all A. thaliana genes and determine which categories they fit in.
At_genes_to_GO <- read.table("At_GO/ATH_to_GO_unique.tsv",
                             header=FALSE, sep="\t", stringsAsFactors = FALSE)

GO_slim_to_At <- list()
for(category in names(go_slim_subcategories)) {
  GO_slim_to_At[[category]] <- character()
}

for(i in 1:nrow(At_genes_to_GO)) {

  At_gene <- At_genes_to_GO[i, 1]
  GO_term <- At_genes_to_GO[i, 2]
  
  for(category in names(go_slim_subcategories)) {
   if(GO_term %in% go_slim_subcategories[[category]]$child_go_id) {
     GO_slim_to_At[[category]] <- c(GO_slim_to_At[[category]], At_gene)
   }
  }
}

for(category in names(go_slim_subcategories)) {
  GO_slim_to_At[[category]] <- table(GO_slim_to_At[[category]])
}

# Keep single biological process per A. thaliana gene. If a gene is involved in multiple high-level processes then only keep the one in the smallest category.

category_lengths <- sort(sapply(GO_slim_to_At, length))

At_genes <- At_genes_to_GO$V1[-which(duplicated(At_genes_to_GO$V1))]

At_gene_single_category <- list()

for(category in names(category_lengths)) {
  At_gene_single_category[[category]] <- character()
}

At_gene_single_category[["Unknown biological processes"]] <- character()

for(At in At_genes) {
  
  At_categorized <- FALSE
  
  for(category in names(category_lengths)) {
     if(At %in% names(GO_slim_to_At[[category]])) {
         At_gene_single_category[[category]] <- c(At_gene_single_category[[category]], At)
         At_categorized <- TRUE
         break
     }
  }
  
  if(! At_categorized) {
    At_gene_single_category[["Unknown biological processes"]] <- c(At_gene_single_category[["Unknown biological processes"]], At)
  }
}


# Save output files.
saveRDS(object = GO_slim_to_At, file = "At_GO_RDS/GO_slim_to_At.rds")
saveRDS(object = At_gene_single_category, file = "At_GO_RDS/At_gene_single_category.rds")

