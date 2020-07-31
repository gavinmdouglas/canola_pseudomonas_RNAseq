# Set working directory.
setwd("Dropbox/work/Langille/pseudomonas/canola_pseudomonas/")

# Read in blast top hits.
blast_top_hits <- read.table("tables/blastn_out_napus_vs_thaliana_evalue0.0001_clean_tophits.txt",
                             header=T, stringsAsFactors = F, sep="\t")

# Remove blast columns that wont be merged.
blast_top_hits <- blast_top_hits[, -c(4,5,6,7,8,9,10)]

# Read in tables of info parsed from B. napus and A. thaliana cds fastas.
napus_info <- read.table("tables/Brassica_napus.AST_PRJEB5043_v1.cds.info.txt",
                         header=T, stringsAsFactors = F, sep="\t", comment.char="",
                         na.strings=c(""), quote="")

thaliana_info <- read.table("tables/Arabidopsis_thaliana.TAIR10.cds.info.txt",
                         header=T, stringsAsFactors = F, sep="\t", comment.char="",
                         na.strings=c(""), quote="")

# Add species name in front of each species annotation column so they can be distinguished.
colnames(napus_info) <- paste("Bnapus", colnames(napus_info), sep = "_") 
colnames(thaliana_info) <- paste("Athaliana", colnames(thaliana_info), sep = "_") 

# Read in past functional annotations (GO and INTERPRO domains).
go_info <- read.table("B.napus_prior_annotation/Brassica_napus_GO", header=F, 
                      sep=" ", stringsAsFactors = F)

ipr_info <- read.table("B.napus_prior_annotation/Brassica_napus_IPR.withdescription", 
                       header=F, sep="\t", stringsAsFactors = F, quote="")

# Set column names for these DFs.
colnames(go_info) <- c("Bnapus_gene_id", "GO_id")
colnames(ipr_info) <- c("Bnapus_gene_id", "INTERPRO_id", "INTERPRO_description")

# Each gene / function combination is on separate lines currently. 
# They can be combined into a single row per gene with "aggregate".
go_info_cat <- aggregate(GO_id ~ Bnapus_gene_id, go_info, paste, collapse = ";")

# Note ipr descriptions weren't included since they were going to add a lot of text.
ipr_info_id_cat <- aggregate(INTERPRO_id ~ Bnapus_gene_id, ipr_info, paste, collapse = ";")


# These two tables can be merged together:
go_and_ipr_info <- merge(go_info_cat, 
                         ipr_info_id_cat, 
                         by="Bnapus_gene_id",
                         all.x=TRUE, 
                         all.y=TRUE)


# Read in diff in top 2 bitscores (will be NA if not applicable).
bitscore_diff <- read.table("tables/blastn_out_napus_vs_thaliana_top2_bitscore_diffs.txt",
                            header=T, stringsAsFactors = F, sep="\t") 


# Merge all of the above tables together into one descriptive table for all Brassica napus genes.
# One table will be added at a time, starting with the napus info df and the blast top hits table.
napus_info_w_blast_top_hits <- merge(napus_info, 
                                     blast_top_hits, 
                                     by.x = "Bnapus_gene_name", 
                                     by.y = "qseqid",
                                     all.x = TRUE,
                                     all.y = TRUE)

# Then merge in diff in bitscore for each Bnapus gene.
napus_info_w_blast_top_hits_w_bitscore_diffs <- merge(napus_info_w_blast_top_hits, 
                                                      bitscore_diff, 
                                                      by.x = "Bnapus_gene_name", 
                                                      by.y = "Bnapus_genes",
                                                      all.x = TRUE,
                                                      all.y = TRUE)

# Then merge in A thaliana annotation info, but do not retain info for genes that weren't blasted.
napus_info_w_blast_top_hits_w_bitscore_diffs_w_thaliana_info <- merge(napus_info_w_blast_top_hits_w_bitscore_diffs, 
                                                                      thaliana_info, 
                                                                      by.x = "sseqid", 
                                                                      by.y = "Athaliana_gene_name",
                                                                      all.x = TRUE,
                                                                      all.y = FALSE)

# Finally merge in GO and IPR annotation ids for Bnapus genes. Again only for gene names that are already in table.
combined_func_table <- merge(napus_info_w_blast_top_hits_w_bitscore_diffs_w_thaliana_info, 
                             go_and_ipr_info, 
                             by.x = "Bnapus_gene_symbol", 
                             by.y = "Bnapus_gene_id",
                             all.x = TRUE,
                             all.y = FALSE)

# Can remove unnecessary columns and reorder all columns.
col_order <- c("Bnapus_gene_name", "Bnapus_gene_symbol", "Bnapus_gene", 
               "Bnapus_position", "Bnapus_description", "Athaliana_gene", 
               "Athaliana_gene_symbol", "Athaliana_position", 
               "Athaliana_description", "pident", "evalue", "bitscore", 
               "top2_bitscore_diff", "GO_id", "INTERPRO_id")

combined_func_table <- combined_func_table[,col_order]

# Output this combined table of all annotation.
write.table(x = combined_func_table, file = "tables/Bnapus_merged_func_annot.txt", 
            col.names=TRUE, row.names=FALSE, quote = FALSE, sep="\t")