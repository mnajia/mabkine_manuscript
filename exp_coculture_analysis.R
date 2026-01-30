#Project: Romee Lab Mabkine Manuscript
#Author: Mohamad Najia
#Objective: analyze RNA-seq from exp1: mabkine stimulation of NK cells cocultured with SKBR3 cells

library(DESeq2)
library(IHW)
library(jsonlite)
library(dplyr)
library(data.table)
library(tximport)
library(rhdf5)
library(ggplot2)
library(BuenColors)
library(ggrepel)
library(scales)
library(rjson)
library(ComplexHeatmap)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(irlba)
library(matrixStats)
library(M3C)
library(clusterProfiler)
library(enrichplot)



####################################################
#Function Declarations
####################################################

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}


####################################################
#Set-up Environment
####################################################

#declare colormaps
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', 
                               '#88CEEF', '#C1D5DC', '#EAD397', 
                               '#FDB31A','#E42A2A', '#A31D1D'))(100)

#initialize variables
project_dir <- "/Volumes/blainey_lab/Mo/exps/20251010_MN_ACS_NZ_SS_novaseq/mabkine_manuscript/"
output_dir <- paste0(project_dir, "rna_seq_analysis/exp_coculture/")
tsv_dir <- paste0(project_dir, "kallisto_output_exp_coculture/")

#set up samplesheet
dirs <- list.dirs(tsv_dir, recursive = FALSE, full.names = FALSE)
tmp <- str_split_fixed(dirs, pattern = "_", n=4)
df.samples <- data.frame(sample_name = dirs,
                         timepoint = tmp[,2],
                         donor = tmp[,3],
                         treatment = tmp[,4])
rownames(df.samples) <- df.samples$sample_name

#import hg19 USCS genomeStudio gene symbol to transcript id mappings
fn <- paste0(project_dir, "/kallisto_custom_index/hg19.annot.cdna.HER2.CAR.gene.symbol.transcript.map")
g2t_map <- fread(fn, data.table = FALSE)
colnames(g2t_map) <- c("transcript_id", "gene_symbol")

#import kallisto abundances
tsv_files <- paste0(tsv_dir, df.samples$sample_name, "/KALLISTO/abundance.tsv")
names(tsv_files) <- df.samples$sample_name
txi <- tximport(tsv_files, type = "kallisto", tx2gene = g2t_map)

#export kallisto matrices 
write.table(txi$abundance,
            file = paste0(output_dir, "kallisto_rna_seq_TPM_matrix_exp_coculture.tsv"), 
            row.names = TRUE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

write.table(txi$counts,
            file = paste0(output_dir, "kallisto_rna_seq_counts_matrix_exp_coculture.tsv"), 
            row.names = TRUE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)



####################################################
#DEG Analysis
####################################################

#set up DESeq 
sampleTable <- df.samples[,c("donor", "treatment"), drop = FALSE]

dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ donor + treatment)
dds <- DESeq(dds)

#perform PCA
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds.keep <- dds[keep,]
nrow(dds.keep)

rld <- rlog(dds.keep, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = c("treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$donor <- sampleTable$donor

pdf(file = paste(output_dir, "exp_stim_pca_rlog_trandformed_RNA_seq_data.pdf", sep = "/"), width = 4.5, height = 3, useDingbats = FALSE)  
ggplot(pcaData, aes(x=PC1, y=PC2, fill=treatment, color=treatment)) + 
  geom_point(aes(shape = factor(donor))) +
  xlab(paste0("PC1: (", percentVar[1], "% Variance Explained)")) +
  ylab(paste0("PC2: (", percentVar[2], "% Variance Explained)")) +
  ggtitle("Exp2 PCA") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())
dev.off()

#extract comparisons for differential expression analysis
comp <- data.frame(Var1 = "Mabkine",
                   Var2 = c("IL-15-1nM", "IL-15-1ng-mL", "IL-15-Monalizumab", "Monalizumab"),
                   cond = "treatment")

for (x in 1:dim(comp)[1]) {
  print(x)
  #extract the DESeq contrast
  resSig <- data.frame(results(dds, contrast=c(comp[x,3], comp[x,1], comp[x,2]), parallel = TRUE, pAdjustMethod = "fdr", alpha = 0.05))
  resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) %>% filter(padj < 0.05) -> out
  out$padj_nlog10 <- -log10(out$padj)
  out$gene <- rownames(out)
  
  #round
  out$baseMean <- round(out$baseMean, 2)
  out$log2FoldChange <- round(out$log2FoldChange, 2)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
  
  #output the differential expression table for statistically significant genes
  write.table(out, 
              file = paste0(output_dir, "DESeq_", comp[x,1], "_v_", comp[x,2], "_Padj05.tsv"), 
              row.names = FALSE, 
              col.names = TRUE, 
              sep = "\t", 
              quote = FALSE)
  
  #output the differential expression table for all genes (significant or not)
  resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) -> out.complete
  out.complete$padj_nlog10 <- -log10(out.complete$padj)
  out.complete$gene <- rownames(out.complete)
  
  write.table(out.complete, 
              file = paste0(output_dir, "DESeq_", comp[x,1], "_v_", comp[x,2], ".tsv"), 
              row.names = FALSE, 
              col.names = TRUE, 
              sep = "\t", 
              quote = FALSE)
  
}


#plot volcanos for each comparison
for (i in 1:dim(comp)[1]) {
  treatment <- comp[i,2]
  fn <- paste0(output_dir, "DESeq_Mabkine_v_", treatment, ".tsv")
  df.plot <- fread(fn, header = TRUE, data.table = FALSE)
  
  df.plot$significant <- FALSE
  df.plot[(df.plot$padj < 0.05 & abs(df.plot$log2FoldChange)>2), "significant"] <- TRUE
  df.plot$label <- ""
  df.plot[df.plot$significant, "label"] <- df.plot[df.plot$significant, "gene"]
  
  gg <- ggplot(df.plot, aes(x=log2FoldChange, y=padj_nlog10, color=significant, alpha=significant)) + 
    geom_point() + 
    scale_color_manual(values=c("#999999", "#E69F00")) + 
    scale_alpha_manual(values=c(0.2,1)) + 
    scale_size_manual(values=c(0.75,1)) + 
    geom_hline(yintercept=(-1*log10(0.05)), linetype="dashed", color = "black") + 
    geom_vline(xintercept=0, linetype="dashed", color = "black") + 
    scale_x_continuous(limits = c(-10, 10)) + 
    scale_y_continuous(limits = c(0, 30)) + 
    xlab("log2(fold change)") + 
    ylab("-log10(P-adjusted)") + 
    ggtitle(paste0("Mabkine v. ", treatment, " DEGs") ) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour="black"),
          legend.position = "none",
          axis.text.y = element_text(colour="black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    geom_text_repel(data = df.plot,
                    aes(label = label), 
                    segment.color = "black",
                    colour = "black",
                    size = 4,
                    box.padding = 1,
                    max.overlaps = 500,
                    show.legend = FALSE)
  
  pdf(file = paste0(output_dir, "exp_stim_DEG_volcano_Mabkine_v_", treatment, ".pdf"), width = 6, height = 6, useDingbats = FALSE)
  print(gg)
  dev.off()
  
}



#plot expression of activation marker genes
df.tpms <- log2(txi$abundance + 1)

activating_genes <- c("CCL3", "CCL4", "CCL5", "CD38", "CSF2", "FASLG", "GZMB", "GZMK", "IFNG", "IL2RA", "LAMP1", "LAMP2", "LTA", "TNF")
df.degs <- df.tpms[activating_genes, c(1,6,11,  5,10,15,  2,7,12,  3,8,13,  4,9,14)]
df.degs.z <- scaleRows(df.degs) %>% as.data.frame()

df.degs.z.avg <- df.degs.z
df.degs.z.avg$IL_15_1ng_mL <- rowMeans(df.degs.z[,c(1:3)])
df.degs.z.avg$Monalizumab <- rowMeans(df.degs.z[,c(4:6)])
df.degs.z.avg$IL_15_1nM <- rowMeans(df.degs.z[,c(7:9)])
df.degs.z.avg$IL_15_Monalizumab <- rowMeans(df.degs.z[,c(10:12)])
df.degs.z.avg$Mabkine <- rowMeans(df.degs.z[,c(13:15)])
df.degs.z.avg[,c(1:15)] <- NULL

ha = rowAnnotation(foo = anno_mark(at = match(activating_genes, rownames(df.degs.z.avg)), labels = activating_genes))

ht1 <- Heatmap(df.degs.z.avg, #df.degs.z
               col = pal_rna, 
               #top_annotation = column_ha, 
               show_row_names = FALSE, 
               show_column_names = TRUE, 
               cluster_rows = TRUE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               clustering_distance_rows = "euclidean", 
               clustering_method_rows = "ward.D2", 
               #row_km = 2,
               row_title = NULL,
               column_title = "NK Activation Marker Genes",
               #row_names_gp = gpar(fontsize = 10),
               #column_names_gp = gpar(fontsize = 10),
               #row_names_side = "left", 
               column_names_side = "top",
               #row_dend_side = "right",
               #rect_gp = gpar(col = "black", lwd = 0.2),
               right_annotation = ha,
               #width = unit(0.5, "cm"),
               heatmap_legend_param = list(direction = "vertical", border = TRUE),
               name = "Z-score of log2(TPM+1)", 
               border = TRUE)

pdf(paste0(output_dir, "exp_coculture_NK_activation_marker_genes_heatmap_donor_avg.pdf"), width = 4, height = 5, useDingbats = FALSE)
draw(ht1, merge_legend = TRUE)
dev.off()



#Update the heatmaps
comp <- data.frame(Var1 = "Mabkine",
                   Var2 = c("IL-15-1nM", "IL-15-1ng-mL", "IL-15-Monalizumab", "Monalizumab"),
                   cond = "treatment")

fn <- paste0(output_dir, "DESeq_", comp[,1], "_v_", comp[,2], "_Padj05.tsv")
res <- lapply(fn, function(x) {
  df <- fread(x, header = TRUE, data.table = FALSE)
  df <- filter(df, abs(log2FoldChange) > 2)
  
  return(df$gene)
})

degs <- unique(unlist(res))

df.tpms <- log2(txi$abundance + 1)
df.degs <- df.tpms[degs, c(1,6,11,  5,10,15,  2,7,12,  3,8,13,  4,9,14)]
df.degs.z <- scaleRows(df.degs) %>% as.data.frame()

df.degs.z.avg <- df.degs.z
df.degs.z.avg$IL_15_1ng_mL <- rowMeans(df.degs.z[,c(1:3)])
df.degs.z.avg$Monalizumab <- rowMeans(df.degs.z[,c(4:6)])
df.degs.z.avg$IL_15_1nM <- rowMeans(df.degs.z[,c(7:9)])
df.degs.z.avg$Monalizumab_IL_15 <- rowMeans(df.degs.z[,c(10:12)])
df.degs.z.avg$Mabkine <- rowMeans(df.degs.z[,c(13:15)])
df.degs.z.avg[,c(1:15)] <- NULL

highlight <- c("ZBTB32","IL13","LTA","IFNG","SLAMF1","CCL3","IL3","CISH","CCL4",
               "IL2RA","GZMB","CCL1","TNF","BATF3","LTB","TNFRSF8",#"IFI30","TEAD4",
               "MYC","TNFRSF12A","DDX21","PRMT1","PRDX1","FOSB","DDX39A",#"TGFBR3L",
               "LAG3","KLRF1","IL7R","ID3","PRDM11","ZBP1","KLRC4","KLRB1","SIGLEC17P",
               "CCR3","ZBTB20","KLF1","ATF3","CX3CR1","CCR4","CCL22","HMOX1",#"DLL3",
               "FGFBP2","TCF7","FCRL3","ADAMTS13","CD27","CD79A","CXCL13") #"BCL6",

highlight <- unique(c(highlight, activating_genes))

ha = rowAnnotation(foo = anno_mark(at = match(highlight, rownames(df.degs.z.avg)), labels = highlight))

ht1 <- Heatmap(df.degs.z.avg, #df.degs.z
               col = pal_rna, 
               #top_annotation = column_ha, 
               show_row_names = FALSE, 
               show_column_names = TRUE, 
               cluster_rows = TRUE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               clustering_distance_rows = "euclidean", 
               clustering_method_rows = "ward.D2", 
               row_km = 2,
               row_title = NULL,
               column_title = "NK cells 24hr co-culture",
               #row_names_gp = gpar(fontsize = 10),
               #column_names_gp = gpar(fontsize = 10),
               #row_names_side = "left", 
               column_names_side = "top",
               #row_dend_side = "right",
               #rect_gp = gpar(col = "black", lwd = 0.2),
               right_annotation = ha,
               #width = unit(0.5, "cm"),
               heatmap_legend_param = list(direction = "vertical", border = TRUE),
               name = "Z-score of log2(TPM+1)", 
               border = TRUE)

pdf(paste0(output_dir, "exp_coculture_DEGs_heatmap_donor_avg.pdf"), width = 5, height = 10, useDingbats = FALSE)
draw(ht1, merge_legend = TRUE)
dev.off()





