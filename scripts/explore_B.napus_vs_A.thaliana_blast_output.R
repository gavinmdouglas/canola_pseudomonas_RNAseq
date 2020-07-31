# This script was used to explore the raw BLASTn output, to output a new table
# with just the top hit per B. napus gene, and to output the difference in bitscores
# between the top 2 A. thaliana genes for cases where multiple genes were hit.

setwd("Dropbox/work/Langille/pseudomonas/canola_pseudomonas/")
blast_out <- read.table("B.napus_vs_A.thaliana_blast_out/blastn_out_napus_vs_thaliana_evalue0.0001.txt", 
                        header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Added column names.
colnames(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                         "gapopen", "qstart", "qend", "sstart", "send", 
                         "evalue", "bitscore")

napus_genes <- unique(blast_out$qseqid)

length(napus_genes)
#[1] 64996

# There were 101040 B. napus gene queried, which means 64.33% of genes had hits.

# Retain only the best blast hit for each AT gene (i.e. get rid of other 
# transcripts).

# Inititalize same dataframe with no columns.
blast_out_subset <- data.frame(matrix(NA, nrow=0,ncol=12))
colnames(blast_out_subset) <- colnames(blast_out)

# Loop over all B. napus genes.
for(gene in napus_genes){
  
  # Get A. thaliana hits for each gene and remove transcript numbers.
  blast_out_gene_subset <- blast_out[which(blast_out$qseqid==gene), , drop=F]
  blast_out_gene_subset$sseqid <- gsub(".\\d+$", "", blast_out_gene_subset$sseqid)
  
  # Get unique hits.
  subset_unique_hits <- unique(blast_out_gene_subset$sseqid)
  
  for (hit in subset_unique_hits){
    blast_out_gene_hit_subset <- blast_out_gene_subset[
                                which(blast_out_gene_subset$sseqid==hit), , drop=F]
  
    max_bit_score <- max(blast_out_gene_hit_subset$bitscore)
    
    blast_out_gene_hit_subset_single <- blast_out_gene_hit_subset[sample(which(blast_out_gene_hit_subset$bitscore == max_bit_score), 1) , , drop=F]
    
    blast_out_subset <- rbind(blast_out_subset, blast_out_gene_hit_subset_single)
    
  }
}

# Plot number of unique A. thaliana gene hits.
num_unique_hits <- c()

for(gene in napus_genes){
  num_unique_hits <- c(num_unique_hits,
                       nrow(blast_out_subset[which(blast_out_subset$qseqid==gene), , drop=F]))
}

hist(num_unique_hits, 
     breaks=1000, 
     xlab="Number of unique A. thaliana genes hit", 
     ylab="Number of B. napus genes", 
     main="",
     xlim=c(0,20),
     ylim=c(0,75000))

# Get distribution of differences in bitscores for top 2 hits for genes with multiple hits.
napus_genes_mult_hit <- napus_genes[which(num_unique_hits > 1)]
blast_out_subset_mult_hit <- blast_out_subset[which(blast_out_subset$qseqid %in% napus_genes_mult_hit), ]

# In same loop get new dataframe with only top hit for each B. napus gene.
napus_genes_single_hit <- napus_genes[which(num_unique_hits == 1)]
blast_out_subset_top_hit <- blast_out_subset[which(blast_out_subset$qseqid %in% napus_genes_single_hit), ]

top2_bitscore_diff <- c()

for(gene in napus_genes_mult_hit){
  
  # Get subset for this gene
  blast_out_subset_mult_hit_gene <- blast_out_subset_mult_hit[which(blast_out_subset_mult_hit$qseqid==gene), ]
  
  # Get ordered bitscores and diff between top 2.
  bitscore_order <- order(blast_out_subset_mult_hit_gene[which(blast_out_subset_mult_hit_gene$qseqid==gene), "bitscore"], decreasing=TRUE)
  mult_hit_top2_bitscores <- blast_out_subset_mult_hit_gene$bitscore[bitscore_order][1:2]
  top2_bitscore_diff <- c(top2_bitscore_diff, mult_hit_top2_bitscores[1]-mult_hit_top2_bitscores[2])
  
  # Based on these ordered bitscores also return the top A. thaliana get hit and add to top hit df.
  blast_out_subset_top_hit <- rbind(blast_out_subset_top_hit, blast_out_subset_mult_hit_gene[bitscore_order[1], ])
  
}

hist(top2_bitscore_diff,
     breaks=100,
     col="grey",
     main="",
     ylab="Number of B. napus genes",
     xlab="Difference in bitscore between top 2 A. thaliana hits",
     xlim=c(0,11000),
     ylim=c(0,5000))
length(which(top2_bitscore_diff < 100))/length(top2_bitscore_diff)
#[1] 0.2407704

# Output cleaned blast output file with only the top transcript shown for each At gene.
write.table(x = blast_out_subset, file = "tables/blastn_out_napus_vs_thaliana_evalue0.0001_clean.txt", 
            col.names=TRUE, row.names=FALSE, quote = FALSE, sep="\t")

# Output only the top hit for each B. napus gene.
write.table(x = blast_out_subset_top_hit, file = "tables/blastn_out_napus_vs_thaliana_evalue0.0001_clean_tophits.txt", 
            col.names=TRUE, row.names=FALSE, quote = FALSE, sep="\t")

# Output the difference in bit scores for each gene with multi-hits.
top2_bitscore_diff_df <- data.frame(matrix(NA,nrow=64996, ncol=2))
colnames(top2_bitscore_diff_df) <- c("Bnapus_genes", "top2_bitscore_diff")
top2_bitscore_diff_df$Bnapus_genes <- c(napus_genes_mult_hit, napus_genes_single_hit)
top2_bitscore_diff_df$top2_bitscore_diff <- c(top2_bitscore_diff, rep(NA, length(napus_genes_single_hit)))

write.table(x = top2_bitscore_diff_df, file = "tables/blastn_out_napus_vs_thaliana_top2_bitscore_diffs.txt", col.names=TRUE, row.names=FALSE, quote = FALSE, sep="\t")
