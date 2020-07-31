library("VennDiagram")

# Function to take in a character vector and will create output file with 1 character per line. 
# Output file name is based on input strings given in separate vector.
write_out_vec <- function(vec2write, str_for_file) {
 
  file_out <- paste(str_for_file, collapse="_")
  print(paste("Writing:", file_out))
  write.table(vec2write, file = file_out, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

# Function that will draw three way Venn diagram and output gene ids overlapping in each set to textfiles.
deseq2_ThreeWayVenn_and_set <- function(results1, results2, results3, name1, name2, name3, plot_outdir,
                                        padj_cut=0.1, l2fc_cut=0, compare_type="de", 
                                        prefix = "genes", suffix="txt",
                                        venn_col=c("#1f78b4", "#33a02c", "#e31a1c")) {
  
  if(compare_type == "de") {
      # Just looking for differentially expressed. Note takes absolute of log2foldchanges.
      set1 <- rownames(results1)[which(results1$padj < padj_cut & abs(results1$log2FoldChange) > l2fc_cut)]
      set2 <- rownames(results2)[which(results2$padj < padj_cut & abs(results2$log2FoldChange) > l2fc_cut)]
      set3 <- rownames(results3)[which(results3$padj < padj_cut & abs(results3$log2FoldChange) > l2fc_cut)]
  } else if(compare_type == "up") {
    # Get up-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange > l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange > l2fc_cut)]
    set3 <- rownames(results3)[which(results3$padj < padj_cut & results3$log2FoldChange > l2fc_cut)]
 
  } else if(compare_type == "down") {
    # Get down-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange < l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange < l2fc_cut)]
    set3 <- rownames(results3)[which(results3$padj < padj_cut & results3$log2FoldChange < l2fc_cut)]
  } else {
    stop("Option compare_type needs to be one of \"de\", \"up\", or \"down\".")
  }
  
  # Get ids overlapping between each set.
  set_n12=set1[which(set1 %in% set2)]
  set_n23=set2[which(set2 %in% set3)]
  set_n13=set1[which(set1 %in% set3)]
  set_n123=set_n12[which(set_n12 %in% set3)]
  
  venn.plot <- draw.triple.venn(area1=length(set1),
                                area2=length(set2),
                                area3=length(set3),
                                n12=length(set_n12),
                                n23=length(set_n23),
                                n13=length(set_n13),
                                n123=length(set_n123),
                                category=c(name1, name2, name3),
                                fill=venn_col)

  out_info <- paste(compare_type, "padj", as.character(padj_cut), "l2fc", as.character(l2fc_cut), sep="_")
  out_suffix <- paste(out_info, suffix, sep=".")
  
  out_pdf <- paste(plot_outdir, "/", basename(prefix), "_", out_info, ".pdf", sep="")
  pdf(file=out_pdf)
  grid.draw(venn.plot);
  dev.off()
  grid.newpage();
  
  # Write out sets to files.
  write_out_vec(vec2write = set1, str_for_file = c(prefix, name1, out_suffix))
  write_out_vec(vec2write = set2, str_for_file = c(prefix, name2, out_suffix))
  write_out_vec(vec2write = set3, str_for_file = c(prefix, name3, out_suffix))
  write_out_vec(vec2write = set_n12, str_for_file = c(prefix, name1, name2, out_suffix))
  write_out_vec(vec2write = set_n23, str_for_file = c(prefix, name2, name3, out_suffix))
  write_out_vec(vec2write = set_n13, str_for_file = c(prefix, name1, name3, out_suffix))
  write_out_vec(vec2write = set_n123, str_for_file = c(prefix, name1, name2, name3, out_suffix))  
}

# Function that will draw two way Venn diagram and output gene ids overlapping in each set to textfiles.
deseq2_TwoWayVenn_and_set <- function(results1, results2, name1, name2, plot_outdir,
                                        padj_cut=0.1, l2fc_cut=0, compare_type="de", 
                                        prefix = "genes", suffix="txt",
                                        venn_col=c("#1f78b4", "#33a02c")) {
  
  if(compare_type == "de") {
    # Just looking for differentially expressed. Note takes absolute of log2foldchanges.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & abs(results1$log2FoldChange) > l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & abs(results2$log2FoldChange) > l2fc_cut)]
  } else if(compare_type == "up") {
    # Get up-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange > l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange > l2fc_cut)]
    
  } else if(compare_type == "down") {
    # Get down-regulated only.
    set1 <- rownames(results1)[which(results1$padj < padj_cut & results1$log2FoldChange < l2fc_cut)]
    set2 <- rownames(results2)[which(results2$padj < padj_cut & results2$log2FoldChange < l2fc_cut)]
  } else {
    stop("Option compare_type needs to be one of \"de\", \"up\", or \"down\".")
  }
  
  # Get ids overlapping between each set.
  set_n12=set1[which(set1 %in% set2)]
  
  if(length(set1) < length(set2)) {
    rotation_setting <- 180
    names_set <- c(name2, name1)
  } else {
    rotation_setting <- 0
    names_set <- c(name1, name2)
  }

  venn.plot <- draw.pairwise.venn(
    area1=length(set1),
    area2=length(set2),
    cross.area=length(set_n12),
    category=names_set,
    cex = c(2, 2, 2),
    cat.cex = c(2, 2),
    fill=venn_col,
    rotation.degree=rotation_setting)
  
  out_info <- paste(compare_type, "padj", as.character(padj_cut), "l2fc", as.character(l2fc_cut), sep="_")
  out_suffix <- paste(out_info, suffix, sep=".")
  
  out_pdf <- paste(plot_outdir, "/", basename(prefix), "_", out_info, ".pdf", sep="")
  pdf(file=out_pdf)
  grid.draw(venn.plot);
  dev.off()
  grid.newpage();
  
  # Write out sets to files.
  write_out_vec(vec2write = set1, str_for_file = c(prefix, name1, out_suffix))
  write_out_vec(vec2write = set2, str_for_file = c(prefix, name2, out_suffix))
  write_out_vec(vec2write = set_n12, str_for_file = c(prefix, name1, name2, out_suffix))
}


