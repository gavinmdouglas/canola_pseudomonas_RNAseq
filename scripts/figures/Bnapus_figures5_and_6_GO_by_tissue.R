# Figures of top GO enrichments based on (1) up-regulated and (2) down-regulated genes.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

library(cowplot)
library(ggplot2)
library(Hmisc)
library(stringr)


go_enrichment_barplot <- function(go_output, title, num_to_keep=15, go2ignore=c(), x_max=13) {
  
  # Remove and GO ids that are specified (which would be because they are 100% redundant with another GO in the plot already)

  row2remove <- which(go_output$filt_results$GO.ID %in% go2ignore)
  if(length(row2remove) > 0) {
    go_output$filt_results <- go_output$filt_results[-row2remove, ]
  }
  
  if(nrow(go_output$filt_results) > num_to_keep) {
    go_output$filt_results <- go_output$filt_results[1:num_to_keep, ]
  }
  
  go_output$filt_results$capitalized <- capitalize(go_output$filt_results$Term)
  
  go_output$filt_results$capitalized <- str_wrap(go_output$filt_results$capitalized, 30)
  
  go_output$filt_results$capitalized <- factor(go_output$filt_results$capitalized,
                                               levels=rev(go_output$filt_results$capitalized))
  return(
    ggplot(go_output$filt_results, aes(x = capitalized, y = fold, fill=log10_fdr)) +
      geom_bar(stat="identity") +
      geom_hline(yintercept=1, linetype="dotted") +
      ylim(0, x_max) +
      coord_flip() +
      xlab("GO annotation") +
      ylab("Fold enrichment") +
      labs(fill=expression('log'[10]*'(q)')) +
      scale_fill_gradient(low = "dark red", high = "orange", limits = c(-30, log10(0.05))) +
      ggtitle(title) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.y = element_text(size=10))
  )
}

go_enrichment_results <- readRDS("At_GO_RDS/go_enrichment_results.rds")

# A lot of terms are cut-off in the dataframe so this list is a dictionary to correct them.
terms2replace <- list()
terms2replace[["generation of precursor metabolites and ..."]] <- "Generation of precursor metabolites and energy"
terms2replace[["regulation of jasmonic acid mediated sig..."]] <- "Regulation of jasmonic acid mediated signaling pathway"
terms2replace[["benzene-containing compound metabolic pr..."]] <- "Benzene-containing compound metabolic process"
terms2replace[["defense response, incompatible interacti..."]] <- "Defense response, incompatible interaction"
terms2replace[["photosynthesis, light harvesting in phot..."]] <- "Photosynthesis, light harvesting in photosystem I"
terms2replace[["cellular polysaccharide metabolic proces..."]] <- "Cellular polysaccharide metabolic process"
terms2replace[["plant-type cell wall organization or bio..."]] <- "Plant-type cell wall organization or biogenesis"
terms2replace[["external encapsulating structure organiz..."]] <- "External encapsulating structure organization"
terms2replace[["branched-chain amino acid catabolic proc..."]] <- "Branched-chain amino acid catabolic process"
terms2replace[["glutamine family amino acid catabolic pr..."]] <- "Glutamine family amino acid catabolic process"
terms2replace[["regulation of response to water deprivat..."]] <- "Regulation of response to water deprivation"
terms2replace[["glutamine family amino acid catabolic pr..."]] <- "Glutamine family amino acid catabolic process"
terms2replace[["glucosamine-containing compound cataboli..."]] <- "Glucosamine-containing compound catabolic process"
terms2replace[["aromatic amino acid family catabolic pro..."]] <- "Aromatic amino acid family catabolic process"
terms2replace[["phenol-containing compound metabolic pro..."]] <- "Phenol-containing compound metabolic process"
terms2replace[["energy derivation by oxidation of organi..."]] <- "Energy derivation by oxidation of organic compounds"
terms2replace[["glucosamine-containing compound metaboli..."]] <- "Glucosamine-containing compound metabolic process"
terms2replace[["cell wall macromolecule catabolic proces..."]] <- "Cell wall macromolecule catabolic process"
terms2replace[["defense response by callose deposition i..."]] <- "Defense response by callose deposition in cell wall"
terms2replace[["photosynthetic electron transport in pho..."]] <- "Photosynthetic electron transport in photosystem I"
terms2replace[["regulation of photosynthesis, light reac..."]] <- "Regulation of photosynthesis, light reaction"
terms2replace[["unsaturated fatty acid biosynthetic proc..."]] <- "Unsaturated fatty acid biosynthetic process"
terms2replace[["regulation of generation of precursor me..."]] <- "Reg. of generation of precursor met. and energy"
terms2replace[["cellular response to decreased oxygen le..."]] <- "Cellular response to decreased oxygen levels"
terms2replace[["indole-containing compound catabolic pro..."]] <- "Indole-containing compound catabolic process"
terms2replace[["cellular response to phosphate starvatio..."]] <- "Cellular response to phosphate starvation"
terms2replace[["indole-containing compound metabolic pro..."]] <- "Indole-containing compound metabolic process"
terms2replace[["secondary metabolite biosynthetic proces... "]] <- "Secondary metabolite biosynthetic process"
terms2replace[["defense response to fungus, incompatible..."]] <- "Defense response to fungus, incompatible interaction"
terms2replace[["modification of morphology or physiology..."]] <- "Modification of morphology or physiology of other organism"
terms2replace[["carbohydrate derivative catabolic proces..."]] <- "Carbohydrate derivative catabolic process"
terms2replace[["regulation of anion transmembrane transp..."]] <- "Regulation of anion transmembrane transport"
terms2replace[["reactive oxygen species metabolic proces..."]] <- "Reactive oxygen species metabolic process"
terms2replace[["secondary metabolite biosynthetic proces..."]] <- "Secondary metabolite biosynthetic process"
terms2replace[["branched-chain amino acid biosynthetic p..."]] <- "Branched-chain amino acid biosynthetic process"




for(category in names(go_enrichment_results)) {
  for(str2replace in names(terms2replace)) {
    if(str2replace %in% go_enrichment_results[[category]]$filt_results$Term) {
      matching_row <- which(go_enrichment_results[[category]]$filt_results$Term == str2replace)
      go_enrichment_results[[category]]$filt_results[matching_row, "Term"] <- terms2replace[[str2replace]]
    }
  }
}


GO_to_AT_table <- read.table("At_GO/ATH_to_GO_unique.tsv", header=FALSE, sep="\t", stringsAsFactors = FALSE)

# Plot barplots of top 15 GO categories while excluding redundant categories.
root_day1_up <- go_enrichment_barplot(go_enrichment_results$Root_up_day1,
                                      title="Root Day 1 Up-regulated",
                                      go2ignore = c("GO:0071453", "GO:0046217", "GO:0052314",
                                                    "GO:0052315", "GO:0009310", "GO:0042402", "GO:0052317"),
                                      x_max = 20)

root_day3_up <- go_enrichment_barplot(go_enrichment_results$Root_up_day3,
                                      title="Root Day 3 Up-regulated",
                                      go2ignore = c("GO:1905622", "GO:0046217", "GO:0052314", "GO:0052315",
                                                    "GO:0052317", "GO:0006030", "GO:0006032", "GO:0046348",
                                                    "GO:1901072", "GO:0071453"),
                                      x_max = 20)  

root_day5_up <- go_enrichment_barplot(go_enrichment_results$Root_up_day5,
                                      title="Root Day 5 Up-regulated",
                                      go2ignore = c(),
                                      x_max = 20)                   

shoot_day1_up <- go_enrichment_barplot(go_enrichment_results$Shoot_up_day1,
                                      title="Shoot Day 1 Up-regulated",
                                      go2ignore = c("GO:0006030", "GO:0006032", "GO:0046348",
                                                    "GO:1901072", "GO:0042343", "GO:0019758",
                                                    "GO:0019761"),
                                      x_max = 20)



shoot_day3_up <- go_enrichment_barplot(go_enrichment_results$Shoot_up_day3,
                                      title="Shoot Day 3 Up-regulated",
                                      go2ignore = c("GO:0071453"),
                                      x_max = 20)

shoot_day5_up <- go_enrichment_barplot(go_enrichment_results$Shoot_up_day5,
                                      title="Shoot Day 5 Up-regulated",
                                      go2ignore = c("GO:0006030", "GO:0006032", "GO:0046348",
                                                    "GO:1901072", "GO:0052317", "GO:0006560",
                                                    "GO:0046217", "GO:0052314", "GO:0052315",
                                                    "GO:1901071"),
                                      x_max = 20)


# Downregulated plots
root_day1_down <- go_enrichment_barplot(go_enrichment_results$Root_down_day1,
                                      title="Root Day 1 Down-regulated",
                                      go2ignore = c("GO:0072347", "GO:0006833", "GO:0018119"),
                                      x_max = 25)

root_day3_down <- go_enrichment_barplot(go_enrichment_results$Root_down_day3,
                                      title="Root Day 3 Down-regulated",
                                      go2ignore = c(),
                                      x_max = 25)                                     

root_day5_down <- go_enrichment_barplot(go_enrichment_results$Root_down_day5,
                                      title="Root Day 5 Down-regulated",
                                      go2ignore = c("GO:0006833"),
                                      x_max = 25)

shoot_day1_down <- go_enrichment_barplot(go_enrichment_results$Shoot_down_day1,
                                       title="Shoot Day 1 Down-regulated",
                                       go2ignore = c(),
                                       x_max = 25)

shoot_day3_down <- go_enrichment_barplot(go_enrichment_results$Shoot_down_day3,
                                       title="Shoot Day 3 Down-regulated",
                                       go2ignore = c(),
                                       x_max = 25)                                     

shoot_day5_down <- go_enrichment_barplot(go_enrichment_results$Shoot_down_day5,
                                       title="Shoot Day 5 Down-regulated",
                                       go2ignore = c("GO:0018119", "GO:0043155", "GO:0006833"),
                                       x_max = 25)  

# Write out plots
pdf(file = "plots/main/Figure5_root_GO.pdf", width=18, height=10, onefile=FALSE)
plot_grid(root_day1_up, root_day3_up, root_day5_up, root_day1_down, root_day5_down,
          nrow=2,
          labels=c('A', 'B', 'C', 'D', 'E'))
dev.off()

tiff(file = "plots/main/Figure5_root_GO.tiff", width=18, height=10,
     compression = "lzw", res=300, units="in")
plot_grid(root_day1_up, root_day3_up, root_day5_up, root_day1_down, root_day5_down,
          nrow=2,
          labels=c('A', 'B', 'C', 'D', 'E'))
dev.off()




pdf(file = "plots/main/Figure6_shoot_GO.pdf", width=12, height=10, onefile=FALSE)
plot_grid(shoot_day1_up, shoot_day3_up,
          shoot_day5_up, shoot_day5_down,
          nrow=2, ncol=2,
          labels=c('A', 'B', 'C', 'D'))

dev.off()

tiff(file = "plots/main/Figure6_shoot_GO.tiff", width=12, height=10,
     compression = "lzw", res=300, units="in")
plot_grid(shoot_day1_up, shoot_day3_up,
          shoot_day5_up, shoot_day5_down,
          nrow=2, ncol=2,
          labels=c('A', 'B', 'C', 'D'))
dev.off()