### Commands to run GO enrichment on all gene sets.

rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)
library(topGO)
library(Hmisc)
library(openxlsx)

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

# First get set of background genes and mappings from At orthologs to GO terms.
gene_background <- rownames(read.table("At_deseq2_outfiles/At_homolog_deseq2_log2fold.txt",
                                       header=TRUE, sep="\t", row.names=1, quote="", comment.char="", stringsAsFactors = FALSE))

At_gene_to_GO_raw <- read.delim(file = "At_GO/ATH_to_GO.tsv", header=FALSE, sep="\t", stringsAsFactors = FALSE)

At_gene_to_GO <- list()

for(g in gene_background) {
  At_gene_to_GO[[g]] <- At_gene_to_GO_raw[which(At_gene_to_GO_raw$V1 == g), "V2"]
}

run_topGO_classic_enrich <- function(sig_genes, go_map, background, name, ontology_set="BP", min_sig=5, max_annotated=100000, min_fdr=0.05) {
  
  geneList <- factor(as.integer(background %in% sig_genes))
  names(geneList) <- background
  
  topGOdata <- new('topGOdata',
                   description=name,
                   ontology=ontology_set,
                   allGenes=geneList,
                   nodeSize = 10,
                   annot=annFUN.gene2GO,
                   gene2GO = go_map)
  
  classic_test <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  
  resultFisher <- getSigGroups(topGOdata, classic_test)
  
  allRes <- GenTable(topGOdata, classic = resultFisher, ranksOf = "classic", topNodes=length(score(resultFisher)))
  allRes$fold <- allRes$Significant/allRes$Expected
  allRes$logfold <- log2(allRes$fold)
  allRes$abslogfold <- abs(allRes$logfold)
  allRes$p <- allRes$classic
  allRes$p[which(allRes$p == "< 1e-30")] <- "1e-30"
  allRes$p <- as.numeric(allRes$p)
  allRes$fdr <- p.adjust(allRes$p, "fdr")
  allRes$log10_fdr <- log10(allRes$fdr)

  # The above table would be useful to save for future reference and
  # the below commands will generate the top hits for visualization.
  allRes_abun <- allRes[which(allRes$Significant >= min_sig & allRes$Annotated <= max_annotated & allRes$fdr < min_fdr), ]
  allRes_abun <- allRes_abun[order(allRes_abun$abslogfold, decreasing = TRUE), ]
  
  return(list(topGOdata=topGOdata,
              resultFisher=resultFisher,
              all_results=allRes,
              filt_results=allRes_abun))
}

# Run GO enrichment for every gene set.
go_enrichment_results <- list()

tissues <- c("Root", "Shoot")
directions <- c("up", "down")
days <- c("day1", "day3", "day5")

for(tissue in tissues) {
  for(direction in directions) {
    for(day in days) {
      gene_set_name = paste(tissue, direction, day, sep="_")
      
      if(direction == "up") {
        filename <- paste("sig_gene_sets/", day, "_up/At_", day, "_genes_", tissue, "_up_padj_0.1_l2fc_2.txt", sep="")
      } else if(direction == "down") {
        filename <- paste("sig_gene_sets/", day, "_down/At_", day, "_genes_", tissue, "_down_padj_0.1_l2fc_-2.txt", sep="")
      }
      
      gene_set <- read.table(filename, stringsAsFactors = FALSE)$V1
      
      go_enrichment_results[[gene_set_name]] <- run_topGO_classic_enrich(sig_genes = gene_set,
                                                                         go_map = At_gene_to_GO,
                                                                         background = gene_background,
                                                                         name = gene_set_name)
      
      rm(gene_set_name)
      rm(gene_set)
      rm(filename)

    }
  }
}

saveRDS(object = go_enrichment_results, file="At_GO/go_enrichment_results.rds")

# Write out full GO enrichment as spreadsheet.

full_GO <- list()

for(genes in names(go_enrichment_results)) {
  full_GO[[genes]] <- go_enrichment_results[[genes]]$all_results
}

write.xlsx(full_GO, file = "tables/gene_set_full_GO_enrichment.xlsx")
