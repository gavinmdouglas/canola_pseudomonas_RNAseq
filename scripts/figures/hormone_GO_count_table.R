# Table of defense hormone GO enrichments - for original draft, but possibly better as a figure actually

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

go_enrichment_results <- readRDS("At_GO/go_enrichment_results.rds")


SA_GO <- "GO:0009751"
JA_GO <- "GO:0009753"
ET_GO <- "GO:0009723"
focal_GO <- c(SA_GO, JA_GO, ET_GO)

root_day1_hormones <- go_enrichment_results$Root_up_day1$all_results[which(go_enrichment_results$Root_up_day1$all_results$GO.ID %in% focal_GO), ]
root_day1_hormones$tissue <- "Root"
root_day1_hormones$day <- 1

root_day3_hormones <- go_enrichment_results$Root_up_day3$all_results[which(go_enrichment_results$Root_up_day3$all_results$GO.ID %in% focal_GO), ]
root_day3_hormones$tissue <- "Root"
root_day3_hormones$day <- 3

root_day5_hormones <- go_enrichment_results$Root_up_day5$all_results[which(go_enrichment_results$Root_up_day5$all_results$GO.ID %in% focal_GO), ]
root_day5_hormones$tissue <- "Root"
root_day5_hormones$day <- 5

shoot_day1_hormones <- go_enrichment_results$Shoot_up_day1$all_results[which(go_enrichment_results$Shoot_up_day1$all_results$GO.ID %in% focal_GO), ]
shoot_day1_hormones$tissue <- "Shoot"
shoot_day1_hormones$day <- 1

shoot_day3_hormones <- go_enrichment_results$Shoot_up_day3$all_results[which(go_enrichment_results$Shoot_up_day3$all_results$GO.ID %in% focal_GO), ]
shoot_day3_hormones$tissue <- "Shoot"
shoot_day3_hormones$day <- 3

shoot_day5_hormones <- go_enrichment_results$Shoot_up_day5$all_results[which(go_enrichment_results$Shoot_up_day5$all_results$GO.ID %in% focal_GO), ]
shoot_day5_hormones$tissue <- "Shoot"
shoot_day5_hormones$day <- 5
