#Project: Romee Lab Mabkine Manuscript
#Author: Mohamad Najia
#Objective: analyze RNA-seq from exp1: mabkine stimulation of NK cells

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
output_dir <- paste0(project_dir, "rna_seq_analysis/exp_stim/")
tsv_dir <- paste0(project_dir, "kallisto_output_exp_stim/")

#set up samplesheet
dirs <- list.dirs(tsv_dir, recursive = FALSE, full.names = FALSE)
tmp <- str_split_fixed(dirs, pattern = "_", n=4)
df.samples <- data.frame(sample_name = dirs,
                         timepoint = tmp[,2],
                         donor = tmp[,3],
                         treatment = tmp[,4])
rownames(df.samples) <- df.samples$sample_name

#import hg19 USCS genomeStudio gene symbol to transcript id mappings
fn <- paste0(project_dir, "/kallisto_index/hg19.annot.cdna.gene.symbol.transcript.map")
g2t_map <- fread(fn, data.table = FALSE)
colnames(g2t_map) <- c("transcript_id", "gene_symbol")

#import kallisto abundances
tsv_files <- paste0(tsv_dir, df.samples$sample_name, "/KALLISTO/abundance.tsv")
names(tsv_files) <- df.samples$sample_name
txi <- tximport(tsv_files, type = "kallisto", tx2gene = g2t_map)

#export kallisto matrices 
write.table(txi$abundance,
            file = paste0(output_dir, "kallisto_rna_seq_TPM_matrix_exp_stim.tsv"), 
            row.names = TRUE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

write.table(txi$counts,
            file = paste0(output_dir, "kallisto_rna_seq_counts_matrix_exp_stim.tsv"), 
            row.names = TRUE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)



####################################################
#DEG Analysis
####################################################

#set up DESeq 
sampleTable <- df.samples[,c("donor", "treatment"), drop = FALSE]
sampleTable$treatment <- factor(sampleTable$treatment, levels = c("Mabkine", "IL-15"))

dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ donor + treatment)
dds <- DESeq(dds)

#perform PCA
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds.keep <- dds[keep,]
nrow(dds.keep)

#vsd <- vst(dds.keep, blind = FALSE)
#head(assay(vsd), 3)
#colData(vsd)
#pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)

rld <- rlog(dds.keep, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = c("treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$donor <- sampleTable$donor

pdf(file = paste(output_dir, "exp_stim_pca_rlog_trandformed_RNA_seq_data.pdf", sep = "/"), width = 4.5, height = 3, useDingbats = FALSE)  
ggplot(pcaData, aes(x=PC1, y=PC2, fill=treatment, color=treatment)) + 
  geom_point(aes(shape = factor(donor))) +
  xlab(paste0("PC1: (", percentVar[1], "% Variance Explained)")) +
  ylab(paste0("PC2: (", percentVar[2], "% Variance Explained)")) +
  ggtitle("Exp1 PCA") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())
dev.off()

#extract comparisons for differential expression analysis

#extract the DESeq contrast
resSig <- data.frame(results(dds, contrast=c("treatment", "Mabkine", "IL-15"), parallel = TRUE, pAdjustMethod = "fdr", alpha = 0.05))
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
            file = paste0(output_dir, "DESeq_Mabkine_v_IL-15_Padj05.tsv"), 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

#output the differential expression table for all genes (significant or not)
resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) -> out.complete
out.complete$padj_nlog10 <- -log10(out.complete$padj)
out.complete$gene <- rownames(out.complete)

write.table(out.complete, 
            file = paste0(output_dir, "DESeq_Mabkine_v_IL-15.tsv"), 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)


#plot volcano of DEGs
df.plot <- out.complete
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
  #scale_x_continuous(limits = c(-10, 10)) + 
  #scale_y_continuous(limits = c(0, 30)) + 
  xlab("log2(fold change)") + 
  ylab("-log10(P adjusted)") + 
  ggtitle("Mabkine v. IL-15 DEGs") + 
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

pdf(file = paste0(output_dir, "exp_stim_Mabkine_v_IL-15_DEG_volcano_plot.pdf"), width = 6, height = 6, useDingbats = FALSE)
gg
dev.off()


#plot heatmap of DEGs
degs <- out$gene
df.tpms <- log2(txi$abundance + 1)
df.degs <- df.tpms[degs, ]
df.degs.z <- scaleRows(df.degs) %>% as.data.frame()
df.degs.z <- df.degs.z[,c(1,3,5,2,4,6)]

df.degs.z.avg <- df.degs.z
df.degs.z.avg$IL_15 <- rowMeans(df.degs.z[,c(1:3)])
df.degs.z.avg$Mabkine <- rowMeans(df.degs.z[,c(4:6)])
df.degs.z.avg[,c(1:6)] <- NULL

color_maps <- c("IL-15" = "#76aed0", "Mabkine" = "#a3062a")
column_ha <- HeatmapAnnotation(Treatment = c(rep("IL-15",3), rep("Mabkine",3) ), 
                               col = list(Treatment = color_maps), 
                               annotation_name_side = "left",
                               border = TRUE)

highlight <- c("CCL1","SLAMF1","IL9R","ZBTB32","SETD7","RUNX2","BATF3","TNFSF4",
               "B3GNT5","CD276","HOXA9","ADAMTS10","CD7","PROCR","FAM131A","CXCR2",
               "KLRB1","FCGR3B","SIGLEC17P","FCRL3","GZMM","FCGR3A","FCRL6",
               "SIGLEC9","CD160","CX3CR1","KLF2","CXCR1")

ha = rowAnnotation(foo = anno_mark(at = match(highlight, rownames(df.degs.z)), labels = highlight))

ht1 <- Heatmap(df.degs.z,
               col = pal_rna, 
               top_annotation = column_ha, 
               show_row_names = FALSE, 
               show_column_names = FALSE, 
               cluster_rows = TRUE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               clustering_distance_rows = "euclidean", 
               clustering_method_rows = "ward.D2",
               row_km = 2,
               row_title = NULL,
               column_title = "Day 5 Stimulated NK cells",
               #row_names_gp = gpar(fontsize = 10),
               #column_names_gp = gpar(fontsize = 10),
               #row_names_side = "left", 
               #row_dend_side = "right",
               #rect_gp = gpar(col = "black", lwd = 0.2),
               right_annotation = ha,
               #width = unit(0.5, "cm"),
               heatmap_legend_param = list(direction = "vertical", border = TRUE),
               name = "Z-score of log2(TPM+1)", 
               border = TRUE)

pdf(paste0(output_dir, "exp_stim_DEGs_heatmap.pdf"), width = 5, height = 6, useDingbats = FALSE)
draw(ht1, merge_legend = TRUE)
dev.off()


ht1 <- Heatmap(df.degs.z.avg,
               col = pal_rna, 
               show_row_names = FALSE, 
               show_column_names = TRUE, 
               cluster_rows = TRUE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               clustering_distance_rows = "euclidean", 
               clustering_method_rows = "ward.D2", 
               row_km = 2,
               row_title = NULL,
               column_title = "Day 5 Stimulated NK cells",
               #row_names_gp = gpar(fontsize = 10),
               #column_names_gp = gpar(fontsize = 10),
               column_names_side = "top",
               #row_names_side = "left", 
               #row_dend_side = "right",
               #rect_gp = gpar(col = "black", lwd = 0.2),
               right_annotation = ha,
               #width = unit(0.5, "cm"),
               heatmap_legend_param = list(direction = "vertical", border = TRUE),
               name = "Z-score of log2(TPM+1)", 
               border = TRUE)

pdf(paste0(output_dir, "exp_stim_DEGs_heatmap_donor_avg.pdf"), width = 5, height = 6, useDingbats = FALSE)
draw(ht1, merge_legend = TRUE)
dev.off()


#plot scatterplot of DEGs for Mabkine v IL-15 comparison
df <- out
rownames(df) <- df$gene
df$highlight <- ""
df[df$gene %in% highlight, "highlight"] <- df[df$gene %in% highlight, "gene"]
df$rank <- 1:dim(df)[1]
df$direction <- "Up-regulated"
df[df$log2FoldChange < 0, "direction"] <- "Down-regulated"

gg1 <- ggplot(df, aes(x = rank, y = log2FoldChange, color = direction, label = highlight)) + 
  geom_point(size = 0.75) +
  #pretty_plot(fontsize = 8) + 
  #L_border() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_y_continuous(limits = c(-5,5) ) + 
  scale_color_manual(values = c("#76aed0", "#a3062a")) + 
  xlab("Gene rank") + 
  ylab("log2(fold-change) [Mabkine / IL-15]") + 
  ggtitle("Day 5 DEGs: Mabkine v. IL-15") +
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  geom_text_repel(data = df,
                  aes(label = highlight), 
                  segment.color = "black",
                  colour = "black",
                  size = 4,
                  box.padding = 1,
                  max.overlaps = 500,
                  show.legend = FALSE)

pdf(paste0(output_dir, "exp_stim_DEGs_log2FC_Mabkine_v_IL-15.pdf"), width = 6, height = 6, useDingbats = FALSE)
gg1
dev.off()



#plot expression of inhibitory genes

inhibitory_genes <- c("BTLA", "CD300A", "SIGLEC7", "SIGLEC9", "TIGIT", "KIR2DL1", "KIR2DL3", "KLRD1", "KLRC1", "KLRG1")
df.degs <- df.tpms[inhibitory_genes, c(1,3,5,2,4,6)]
df.degs.z <- scaleRows(df.degs) %>% as.data.frame()

df.degs.z.avg <- df.degs.z
df.degs.z.avg$IL_15 <- rowMeans(df.degs.z[,c(1:3)])
df.degs.z.avg$Mabkine <- rowMeans(df.degs.z[,c(4:6)])
df.degs.z.avg[,c(1:6)] <- NULL


ha = rowAnnotation(foo = anno_mark(at = match(inhibitory_genes, rownames(df.degs.z.avg)), labels = inhibitory_genes))

ht1 <- Heatmap(df.degs.z.avg, #df.degs.z
               col = pal_rna, 
               #top_annotation = column_ha, 
               show_row_names = FALSE, 
               show_column_names = TRUE, 
               cluster_rows = TRUE, 
               cluster_columns = FALSE, 
               show_column_dend = FALSE,
               clustering_distance_rows = "euclidean", #euclidean, pearson
               clustering_method_rows = "ward.D2", #"ward.D2",
               #row_km = 2,
               row_title = NULL,
               column_title = "Inhibitory NK Genes",
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

pdf(paste0(output_dir, "exp_stim_NK_inhibitory_genes_heatmap_donor_avg.pdf"), width = 4, height = 4, useDingbats = FALSE)
draw(ht1, merge_legend = TRUE)
dev.off()



