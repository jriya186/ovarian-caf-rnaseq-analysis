library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggvenn)
library(clusterProfiler)

# Loading counts matrix
counts <- read.table("~/rnaseq-pipeline/results/featurecounts/counts_matrix.txt",
                     header = TRUE,
                     skip = 1,
                     row.names = 1)
dim(counts)
colnames(counts)

# cleaning matrix
count_matrix <- counts[, 6:17]
colnames(count_matrix) <- gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam", "", colnames(count_matrix))
colnames(count_matrix)

# Creating metadata 
metadata <- data.frame(
  sample = colnames(count_matrix),
  condition = c("HUF", "Kuramochi_CAF", "Primary_CAF", 
                "SKOV3_CAF", "SKOV3_CAF", "SKOV3_CAF",
                "Kuramochi_CAF", "Primary_CAF", "Primary_CAF",
                "HUF", "HUF", "Kuramochi_CAF"),
  row.names = colnames(count_matrix)
)

# Setting HUF as reference level
metadata$condition <- factor(metadata$condition,
                             levels = c("HUF","Kuramochi_CAF", "Primary_CAF","SKOV3_CAF"))
table(metadata$condition)

# DESeq2 object
dds<- DESeqDataSetFromMatrix(countData=count_matrix,
                             colData = metadata,
                             design = ~ condition)
dds <- DESeq(dds)

resultsNames(dds)

kuramochi_res <- results(dds,
                         name = "condition_Kuramochi_CAF_vs_HUF",
                         alpha = 0.05)
primary_res <- results(dds,
                       name = "condition_Primary_CAF_vs_HUF",
                       alpha = 0.05)
skov3_res <- results(dds,
                     name = "condition_SKOV3_CAF_vs_HUF",
                     alpha = 0.05)

summary(kuramochi_res)
summary(primary_res)
summary(skov3_res)

# Converting to dataframes
kuramochi_df <- as.data.frame(kuramochi_res)
primary_df <- as.data.frame(primary_res)
skov3_df <- as.data.frame(skov3_res)

# Adding gene symbols to all three
kuramochi_df$symbol <- mapIds(org.Hs.eg.db,
                              keys = rownames(kuramochi_df),
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")

primary_df$symbol <- mapIds(org.Hs.eg.db,
                            keys = rownames(primary_df),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")

skov3_df$symbol <- mapIds(org.Hs.eg.db,
                          keys = rownames(skov3_df),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# Checking top genes for kuramochi
head(kuramochi_df[order(kuramochi_df$padj), c("symbol", "log2FoldChange", "padj")])

EnhancedVolcano(kuramochi_df,
                lab = kuramochi_df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Kuramochi CAF (FC:1) vs HUF',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.3,
                labSize = 3,
                colAlpha = 0.5,
                xlim = c(-8, 8),
                ylim = c(0, 200))

EnhancedVolcano(primary_df,
                lab = primary_df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Primary CAF (FC:1) vs HUF',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.3,
                labSize = 3,
                colAlpha = 0.5,
                xlim = c(-8, 8),
                ylim = c(0, 200))
EnhancedVolcano(skov3_df,
                lab = skov3_df$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'SKOV3 CAF (FC:1) vs HUF',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.3,
                labSize = 3,
                colAlpha = 0.5,
                xlim = c(-8, 8),
                ylim = c(0, 200))

# Significant DEGs for each comparison

kuramochi_sig <- rownames(kuramochi_df[!is.na(kuramochi_df$padj) & 
                                         kuramochi_df$padj < 0.05 & 
                                         abs(kuramochi_df$log2FoldChange) > 1, ])

primary_sig <- rownames(primary_df[!is.na(primary_df$padj) & 
                                     primary_df$padj < 0.05 & 
                                     abs(primary_df$log2FoldChange) > 1, ])

skov3_sig <- rownames(skov3_df[!is.na(skov3_df$padj) & 
                                 skov3_df$padj < 0.05 & 
                                 abs(skov3_df$log2FoldChange) > 1, ])

length(kuramochi_sig)
length(primary_sig)
length(skov3_sig)

# Venn Diagram

deg_list <- list(
  Kuramochi_CAF = kuramochi_sig,
  Primary_CAF = primary_sig,
  SKOV3_CAF = skov3_sig
)

ggvenn(deg_list,
       fill_color = c('#E69F00', '#56B4E9', '#009E73'),
       stroke_size = 0.5,
       set_name_size = 4)

core_caf_genes <- kuramochi_sig

# iterating to keep genes that appear in the next comparison 
for (sig_genes in list(primary_sig, skov3_sig)) {
  core_caf_genes <- core_caf_genes[core_caf_genes %in% sig_genes]
}

length(core_caf_genes)

# adding symbols
core_caf_df <- kuramochi_df[core_caf_genes, c("log2FoldChange", "padj", "symbol")]
core_caf_df <- core_caf_df[order(abs(core_caf_df$log2FoldChange), decreasing = TRUE), ]
head(core_caf_df, 20)

# ranking most significant and most differentially expressed genes
# first replacing 0 padj values with min non-zero padj
min_padj <- min(core_caf_df$padj[core_caf_df$padj > 0], na.rm = TRUE)
core_caf_df$padj_adj <- ifelse(core_caf_df$padj == 0, min_padj, core_caf_df$padj)
# computing scores
core_caf_df$score <- -log10(core_caf_df$padj_adj) * abs(core_caf_df$log2FoldChange)
core_caf_df <- core_caf_df[order(core_caf_df$score, decreasing = TRUE), ]
head(core_caf_df, 20)

# normalization for heatmap
vsd <- vst(dds, blind = FALSE)
vsd_matrix <- assay(vsd)
top50_genes <- rownames(core_caf_df)[!is.na(core_caf_df$symbol[1:55])][1:50]
heatmap_matrix <- vsd_matrix[top50_genes, ]
rownames(heatmap_matrix) <- core_caf_df[top50_genes, "symbol"]
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
# add gene symbols to replace ensmbl ids
rownames(heatmap_matrix) <- core_caf_df[top50_genes, "symbol"]
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
# plot
pheatmap(heatmap_matrix_scaled,
         annotation_col = data.frame(condition = metadata$condition, 
                                     row.names = colnames(heatmap_matrix)),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
         main = "Core CAF Signature — Top 50 Genes")
round(heatmap_matrix_scaled, 2)
## Main Analysis — Resister Genes
## Definition: genes significantly changed in Primary CAF vs HUF (either direction)
## that conditioned CAFs (Kuramochi AND SKOV3) failed to recapitulate

# Step 1: What does full CAF conversion change? (Primary CAF vs HUF)
primary_up <- rownames(primary_df)[
  which(primary_df$padj < 0.05 & primary_df$log2FoldChange > 1)
]
primary_down <- rownames(primary_df)[
  which(primary_df$padj < 0.05 & primary_df$log2FoldChange < -1)
]
primary_changed <- c(primary_up, primary_down)
length(primary_changed)

# Step 2: What did conditioned CAFs change?
kuramochi_changed <- rownames(kuramochi_df)[
  which(kuramochi_df$padj < 0.05 & abs(kuramochi_df$log2FoldChange) > 1)
]
skov3_changed <- rownames(skov3_df)[
  which(skov3_df$padj < 0.05 & abs(skov3_df$log2FoldChange) > 1)
]

# Step 3: Resisters = changed in Primary but missed by BOTH conditioned CAFs
resisters <- primary_changed[
  !primary_changed %in% kuramochi_changed &
    !primary_changed %in% skov3_changed
]
length(resisters)

# Step 4: Inspect resisters
resisters_df <- primary_df[resisters, c("log2FoldChange", "padj", "symbol")]
resisters_df <- resisters_df[order(abs(resisters_df$log2FoldChange), decreasing = TRUE), ]
head(resisters_df, 20)

# Pathway Enrichment

# Convert ENSEMBL IDs to. Entrez IDs

resisters_entrez <- bitr(resisters,
                         fromType = "ENSEMBL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

head(resisters_entrez)
nrow(resisters_entrez)

## Table

# Split resisters into up and down
resisters_up <- resisters_df[resisters_df$log2FoldChange > 0, ]
resisters_down <- resisters_df[resisters_df$log2FoldChange < 0, ]

# Top 10 each, ordered by absolute LFC
top_resisters_up <- head(resisters_up[order(resisters_up$log2FoldChange, decreasing = TRUE), ], 10)
top_resisters_down <- head(resisters_down[order(resisters_down$log2FoldChange), ], 10)

# Add direction column
top_resisters_up$direction <- "Upregulated"
top_resisters_down$direction <- "Downregulated"

# Combine
top_resisters_table <- rbind(top_resisters_up, top_resisters_down)

# Clean up for display
top_resisters_table$log2FoldChange <- round(top_resisters_table$log2FoldChange, 2)
top_resisters_table$padj <- formatC(top_resisters_table$padj, format = "e", digits = 2)

# Display
top_resisters_table[, c("symbol", "log2FoldChange", "padj", "direction")]

# GO Biological Process enrichment

go_resisters <- enrichGO(gene = resisters_entrez$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)

# Checking results
dim(go_resisters@result)
head(go_resisters@result[, c("Description","GeneRatio","pvalue","p.adjust")], 20)

# Barplot
go_results <- go_resisters@result[1:20, ]
go_results$Description <- factor(go_results$Description, 
                                 levels = rev(go_results$Description))

ggplot(go_results, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "firebrick", high = "steelblue") +
  theme_bw() +
  labs(title = "GO Enrichment — Resister Genes",
       x = "Gene Count",
       y = NULL,
       fill = "Adjusted p-value")

# Dotplot
ggplot(go_results, aes(x = Count, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "firebrick", high = "steelblue") +
  scale_size_continuous(range = c(3, 10)) +
  theme_bw() +
  labs(title = "GO Enrichment — Resister Genes",
       x = "Gene Count",
       y = NULL,
       color = "Adjusted p-value",
       size = "Gene Count")

























