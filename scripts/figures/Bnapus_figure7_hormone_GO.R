# Figure of hormone GO enrichments

rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)
library(Hmisc)
library(stringr)
library(ggpubr)

setwd("/home/gavin/projects/pseudomonas_RNAseq/canola_pseudomonas_RNAseq/")

go_enrichment_results <- readRDS("At_GO_RDS/go_enrichment_results.rds")


SA_GO <- "GO:0009751"
JA_GO <- "GO:0009753"
ET_GO <- "GO:0009723"
focal_GO <- c(SA_GO, JA_GO, ET_GO)

root_day1_hormones <- go_enrichment_results$Root_up_day1$all_results[which(go_enrichment_results$Root_up_day1$all_results$GO.ID %in% focal_GO), ]
root_day1_hormones$tissue <- "Root"
root_day1_hormones$Day <- 1

root_day3_hormones <- go_enrichment_results$Root_up_day3$all_results[which(go_enrichment_results$Root_up_day3$all_results$GO.ID %in% focal_GO), ]
root_day3_hormones$tissue <- "Root"
root_day3_hormones$Day <- 3

root_day5_hormones <- go_enrichment_results$Root_up_day5$all_results[which(go_enrichment_results$Root_up_day5$all_results$GO.ID %in% focal_GO), ]
root_day5_hormones$tissue <- "Root"
root_day5_hormones$Day <- 5

shoot_day1_hormones <- go_enrichment_results$Shoot_up_day1$all_results[which(go_enrichment_results$Shoot_up_day1$all_results$GO.ID %in% focal_GO), ]
shoot_day1_hormones$tissue <- "Shoot"
shoot_day1_hormones$Day <- 1

shoot_day3_hormones <- go_enrichment_results$Shoot_up_day3$all_results[which(go_enrichment_results$Shoot_up_day3$all_results$GO.ID %in% focal_GO), ]
shoot_day3_hormones$tissue <- "Shoot"
shoot_day3_hormones$Day <- 3

shoot_day5_hormones <- go_enrichment_results$Shoot_up_day5$all_results[which(go_enrichment_results$Shoot_up_day5$all_results$GO.ID %in% focal_GO), ]
shoot_day5_hormones$tissue <- "Shoot"
shoot_day5_hormones$Day <- 5

root_hormones <- rbind(root_day1_hormones, root_day3_hormones, root_day5_hormones)
root_hormones$Day <- as.factor(root_hormones$Day)
root_hormones$Term <- capitalize(root_hormones$Term)
root_hormones$Term <- str_wrap(root_hormones$Term, 15)

shoot_hormones <- rbind(shoot_day1_hormones, shoot_day3_hormones, shoot_day5_hormones)
shoot_hormones$Day <- as.factor(shoot_hormones$Day)
shoot_hormones$Term <- capitalize(shoot_hormones$Term)
shoot_hormones$Term <- str_wrap(shoot_hormones$Term, 15)

shoot_hormones$p.signif <- "ns"
shoot_hormones[which(shoot_hormones$fdr < 0.1), "p.signif"] <- "*"
shoot_hormones[which(shoot_hormones$fdr < 0.05), "p.signif"] <- "**"
shoot_hormones[which(shoot_hormones$fdr < 0.0001), "p.signif"] <- "**"
shoot_hormones$group1 <- shoot_hormones$Term
shoot_hormones$group2 <- shoot_hormones$Term
shoot_hormones$y.position <- shoot_hormones$fold + 0.2


root_hormones$p.signif <- "ns"
root_hormones[which(root_hormones$fdr < 0.1), "p.signif"] <- "*"
root_hormones[which(root_hormones$fdr < 0.05), "p.signif"] <- "**"
root_hormones[which(root_hormones$fdr < 0.0001), "p.signif"] <- "***"
root_hormones$group1 <- root_hormones$Term
root_hormones$group2 <- root_hormones$Term
root_hormones$y.position <- root_hormones$fold + 0.2

root_hormone_barplot <- ggplot(root_hormones, aes(x = Term, y = fold, fill=Day)) +
                              geom_bar(position = "dodge", stat="identity", colour="black") +
                              geom_hline(yintercept=1, linetype="dotted") +
                              scale_fill_manual(values=c("white", "black", "grey")) +
                              xlab("GO annotation") +
                              ylab("Fold enrichment") +
                              ggtitle("Root") +
                              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                    axis.text.y = element_text(size=10), legend.position = c(0.05, 0.75),
                                    legend.background = element_blank(),
                                    legend.box.background = element_rect(colour = "black")) +
                              scale_y_continuous(expand = c(0, 0), limits = c(0, 7)) +
                              stat_pvalue_manual(data = root_hormones,
                                                 label = "p.signif",
                                                 label.size = 4,
                                                 xmax = NULL,
                                                 position = position_dodge(.90))


shoot_hormone_barplot <- ggplot(shoot_hormones, aes(x = Term, y = fold, fill=Day)) +
                                  geom_bar(position = "dodge", stat="identity", colour="black") +
                                  geom_hline(yintercept=1, linetype="dotted") +
                                  scale_fill_manual(values=c("white", "black", "grey")) +
                                  xlab("GO annotation") +
                                  ylab("Fold enrichment") +
                                  ggtitle("Shoot") +
                                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                        axis.text.y = element_text(size=10), legend.position = c(0.05, 0.75),
                                        legend.background = element_blank(),
                                        legend.box.background = element_rect(colour = "black")) +
                                  scale_y_continuous(expand = c(0, 0), limits = c(0, 7)) +
                                  stat_pvalue_manual(data = shoot_hormones,
                                                     label = "p.signif",
                                                     label.size = 4,
                                                     xmax = NULL,
                                                     position = position_dodge(.90))


pdf(file = "plots/main/Figure7_hormone_GO.pdf", width=12, height=5, onefile=FALSE)
plot_grid(root_hormone_barplot, shoot_hormone_barplot, labels=c('A', 'B'))
dev.off()


tiff(file = "plots/main/Figure7_hormone_GO.tiff",  width=12, height=5,
     compression = "lzw", res=300, units="in")
plot_grid(root_hormone_barplot, shoot_hormone_barplot, labels=c('A', 'B'))
dev.off()
