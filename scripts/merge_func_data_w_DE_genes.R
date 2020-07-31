# Set working directory.
setwd("Dropbox/work/Langille/pseudomonas/canola_pseudomonas/")

# Define function that will take in a filename, a dataframe,the column
# names to merge on, and the output file. Note that the output rows will
# only correspond to the overlapping ids found in the file, not the df.

merge_file_with_df <- function(file_x, df_y, colx, coly, outfile){
  
  in_file <- read.table(file_x, sep="\t", stringsAsFactors = FALSE, 
                      comment.char = "", quote="", header=TRUE)
  
  merged_df <- merge(in_file,
                     df_y,
                     by.x = colx,
                     by.y = coly,
                     all.x = TRUE,
                     all.y = FALSE)
  
  write.table(x = merged_df, file = outfile, 
              col.names=TRUE, row.names=FALSE, quote = FALSE, sep="\t")
}
  
# Read in table of all combined annotations.
combined_annotation <- read.table("tables/Bnapus_merged_func_annot.txt",
                            header=T, stringsAsFactors = F, sep="\t", 
                            comment.char="", quote="")

# Run merge_file_with_df on all original log fold files.
merge_file_with_df(file_x = "RNAseq_DE_genes/original_output/RootDay1.txt",
                   df_y = combined_annotation,
                   colx = "Gene.ID",
                   coly = "Bnapus_gene",
                   outfile = "RNAseq_DE_genes/annotated_output/RootDay1_annot.txt")

merge_file_with_df(file_x = "RNAseq_DE_genes/original_output/RootDay3.txt",
                   df_y = combined_annotation,
                   colx = "Gene.ID",
                   coly = "Bnapus_gene",
                   outfile = "RNAseq_DE_genes/annotated_output/RootDay3_annot.txt")


merge_file_with_df(file_x = "RNAseq_DE_genes/original_output/RootDay5.txt",
                   df_y = combined_annotation,
                   colx = "Gene.ID",
                   coly = "Bnapus_gene",
                   outfile = "RNAseq_DE_genes/annotated_output/RootDay5_annot.txt")


merge_file_with_df(file_x = "RNAseq_DE_genes/original_output/ShootDay1.txt",
                   df_y = combined_annotation,
                   colx = "Gene.ID",
                   coly = "Bnapus_gene",
                   outfile = "RNAseq_DE_genes/annotated_output/ShootDay1_annot.txt")

merge_file_with_df(file_x = "RNAseq_DE_genes/original_output/ShootDay3.txt",
                   df_y = combined_annotation,
                   colx = "Gene.ID",
                   coly = "Bnapus_gene",
                   outfile = "RNAseq_DE_genes/annotated_output/ShootDay3_annot.txt")


merge_file_with_df(file_x = "RNAseq_DE_genes/original_output/ShootDay5.txt",
                   df_y = combined_annotation,
                   colx = "Gene.ID",
                   coly = "Bnapus_gene",
                   outfile = "RNAseq_DE_genes/annotated_output/ShootDay5_annot.txt")