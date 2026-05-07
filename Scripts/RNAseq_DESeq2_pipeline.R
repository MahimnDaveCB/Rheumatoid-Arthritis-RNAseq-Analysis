setwd('add path according to your file organisation')

# Load required libraries
library(DESeq2) # Differential expression analysis
library(pheatmap)  # Heatmap visualization
library(RColorBrewer)  # Color palettes for plots

# Import count matrix
# raw_counts.tsv:
# Rows = genes
# Columns = samples
counts_table <- read.csv('raw_counts.csv', row.names = 1) 
head(counts_table)
dim(counts_table)


# Import sample metadata
# design.tsv should contain:
# sample names + experimental groups
sample_info <- read.csv('design.csv', row.names=1)

# View metadata
sample_info
dim(sample_info)

dim(counts_table)

# Create experimental groups
factors <- factor(sample_info$Group)

# Extract unique groups
groups <- unique(sample_info$Group)

# Reverse group order
groups <- rev(groups)



# IMPORTANT:
# Count matrix IDs and metadata IDs are different
# but sequentially correspond to each other
# =========================================================

# Save original GSM IDs
sample_info$GSM_ID <- rownames(sample_info)

# Replace metadata rownames with counts matrix column names

rownames(sample_info) <- colnames(counts_table)

# Verify matching
all(colnames(counts_table) == rownames(sample_info))

# Expected output:
# TRUE
# Create DESeq2 dataset
# -------------------------
# Design formula specifies
# comparison based on Group
dds <- DESeqDataSetFromMatrix(countData = counts_table, colData=sample_info, design = ~Group)

# Set reference/control group
# -------------------------
# Change "control" if needed
dds$Group <- relevel(dds$Group, ref="undifferentiated_M0_Vehicle")

# Filter low-expression genes
# -------------------------
# Keep genes with >=10 counts in at least minimum group size
keep <- rowSums(counts(dds) >=10) >= min(table(sample_info$Group))

dds <- dds[keep,]

# Run DESeq2 analysis
dds <- DESeq(dds, test="Wald", sfType='poscount')

# Extract differential expression results
deseq_result <- results(dds)

# Convert to dataframe
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)
head(deseq_result)

dim(deseq_result)

names(deseq_result)

# Add gene names as a column
deseq_result$GeneName <- row.names(deseq_result)
names(deseq_result)
head(deseq_result)

# Reorder/select important columns
deseq_result <- subset(deseq_result,
                       select = c("GeneName","padj","pvalue","lfcSE","stat","log2FoldChange","baseMean")
)

names(deseq_result)


# Save complete DESeq2 results
write.table(deseq_result, file='deseq.result.all.tsv', row.names=F, sep='\t')


# Filter significant DEGs
# -------------------------
# padj < 0.05
# |log2FC| >= 1
deg <- subset(deseq_result, padj<0.05 & abs(log2FoldChange)>=1)

dim(deg)  # Number of DEGs
dim(deseq_result) # Total genes analyzed

deg <- deg[order(deg$padj),] # Sort DEGs by adjusted p-value
head(deg)

# Save significant DEGs
write.table(deg, file="deseq_deg.tsv", row.names=F, sep='\t')

# Quality Control & Visualization
# =========================================================

# -------------------------
# Dispersion plot
# -------------------------
# Shows gene-wise dispersion estimates
plotDispEsts(dds, main='GSE328878 Dispersion Estimates')

# Histogram of adjusted p-values
hist(deseq_result$padj, breaks=seq(0,1,length=21), col = 'grey', border = 'white', 
     xlab="", ylab="", ylim=c(0,8000), main='GSE203159 Frequencies of padj-values')

# Volcano Plot
# Define custom colors
old.pal <- palette(c("#00BFFF", "#FF3030"))
# Adjust plot margins
par(mar=c(4,4,2,1), cex.main=1.5)
# Plot title
title=paste(groups[1],"vs",groups[2])
# Base volcano plot
plot(deseq_result$log2FoldChange, -log10(deseq_result$padj), main=title,
     xlab="log2FC",
     ylab="-log10(padj)", pch=20, cex=0.5)

# Highlight significant genes
with( subset( deseq_result, padj <0.05 & abs(log2FoldChange) >=1),
      points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) +3)/2, cex=1))
# Add legend
legend("bottomleft", title=paste("padj<", 0.05, sep=""),
       legend=c("down","up"), pch=20, col=1:2)



# PCA & Heatmap
# =========================================================

# -------------------------
# Variance stabilizing transformation
# -------------------------
# Reduces heteroscedasticity
vsd <- vst(dds,blind=FALSE)

# PCA plot
# -------------------------
# Visualizes sample clustering
plotPCA(vsd, intgroup=c("Group"))

# Extract normalized counts
normalized_counts <- counts(dds, normalized=T)
head(normalized_counts)

# Log2 transformation
# -------------------------
transform_counts <- log2(normalized_counts+1)
head(transform_counts)


# Select top DEGs for heatmap
top_hits <- row.names(deg[1:10,])
head(top_hits)
top_hits

# Extract expression values
top_hits <- transform_counts[top_hits,]
head(top_hits)
#----------------------
# Select top 20 genes
top_genes <- head(deg$Gene, 20)

# Extract transformed expression matrix
mat <- assay(vsd)[top_genes, ]

library(org.Hs.eg.db)
library(AnnotationDbi)

rownames(mat) <- gsub("\\..*", "", rownames(mat))
gene_names <- mapIds(
  org.Hs.eg.db,
  keys = rownames(mat),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
gene_names[is.na(gene_names)] <- rownames(mat)[is.na(gene_names)]
rownames(mat) <- gene_names
colnames(mat) <- paste(
  sample_info$Group,
  sample_info$Donor,
  sep = "_"
)
library(pheatmap)

pheatmap(
  mat,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Top Differentially Expressed Genes"
)

## Heatmap of top DEGs
# -------------------------
#pheatmap(top_hits,cluster_rows = FALSE, cluster_cols=FALSE)


# =========================
# Functional Enrichment Analysis (GO + KEGG)
# =========================


# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)   # Change if organism is not human
library(enrichplot)

# =========================
# Step 1: Prepare gene list
# =========================

gene_symbols <- deg$GeneName

# Convert SYMBOL → ENTREZ ID
gene_df <- bitr(gene_symbols,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# =========================
# Step 2: Background (recommended)
# =========================

universe <- bitr(deseq_result$GeneName,
                 fromType = "ENSEMBL",
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

# =========================
# Step 3: GO Enrichment
# =========================


ego <- enrichGO(gene          = gene_df$ENTREZID,
                universe      = universe$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",   # BP / MF / CC
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# =========================
# Step 4: KEGG Enrichment
# =========================

ekegg <- enrichKEGG(gene         = gene_df$ENTREZID,
                    organism     = "hsa",   # human
                    pvalueCutoff = 0.05)

# =========================
# Step 5: Visualization
# =========================

# GO plots
dotplot(ego, showCategory = 15, title = "GO Biological Process")
barplot(ego, showCategory = 15, title = "GO Enrichment")

# KEGG plots
dotplot(ekegg, showCategory = 15, title = "KEGG Pathways")

# =========================
# Step 6: Save results
# =========================

write.table(as.data.frame(ego),
            file = "GO_enrichment.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(as.data.frame(ekegg),
            file = "KEGG_enrichment.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# KEGG network
cnetplot(ekegg,
         showCategory = 8,
         circular = TRUE,
         colorEdge = TRUE)

library(STRINGdb)

# Initialize STRING
string_db <- STRINGdb$new(version="11.5",
                          species=9606,  # human
                          score_threshold=900)

# Map genes
mapped <- string_db$map(data.frame(gene=deg$GeneName),
                        "gene",
                        removeUnmappedRows = TRUE)

# Get interaction network
hits <- mapped$STRING_id
network <- string_db$get_subnetwork(hits)

# Plot network
string_db$plot_network(hits)

heatplot(ego, showCategory = 10)

selected_genes <- ego@result$geneID[1]  # top GO term genes
selected_genes <- unlist(strsplit(selected_genes, "/"))

mapped_subset <- string_db$map(data.frame(gene=selected_genes),
                               "gene",
                               removeUnmappedRows = TRUE)

string_db$plot_network(mapped_subset$STRING_id)

## RA vs Healthy in M1-like LPS Analysis

selected_samples <- colnames(vsd)[
  sample_info$Group == "M1-like_LPS" &
    sample_info$Donor %in% c("HD1", "RA1")
]
top_genes <- head(deg$GeneName, 30)

mat_subset <- assay(vsd)[top_genes, selected_samples]

library(org.Hs.eg.db)
library(AnnotationDbi)

rownames(mat_subset) <- gsub(
  "\\..*",
  "",
  rownames(mat_subset)
)

gene_names <- mapIds(
  org.Hs.eg.db,
  keys = rownames(mat_subset),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_names[is.na(gene_names)] <-
  rownames(mat_subset)[is.na(gene_names)]

rownames(mat_subset) <- gene_names

colnames(mat_subset) <- paste(
  sample_info$Group[
    colnames(vsd) %in% selected_samples
  ],
  sample_info$Donor[
    colnames(vsd) %in% selected_samples
  ],
  sep = "_"
)

library(pheatmap)

pheatmap(
  mat_subset,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "M1-like LPS: Healthy vs RA"
)


# =========================================================
# Create Disease-Specific Groups
# =========================================================

sample_info$Disease_Group <- paste(
  sample_info$Disease,
  sample_info$Group,
  sep = "_"
)

# View groups
table(sample_info$Disease_Group)

# =========================================================
# Create New DESeq2 Dataset
# =========================================================

dds2 <- DESeqDataSetFromMatrix(
  countData = counts_table,
  colData = sample_info,
  design = ~ Disease_Group
)

# =========================================================
# Set Reference Group
# =========================================================

dds2$Disease_Group <- relevel(
  dds2$Disease_Group,
  ref = "healthy_M1-like_LPS"
)

# =========================================================
# Filter Low Count Genes
# =========================================================

keep <- rowSums(counts(dds2) >= 10) >= 2

dds2 <- dds2[keep, ]

# =========================================================
# Run DESeq2
# =========================================================

dds2 <- DESeq(dds2)

# =========================================================
# View Available Comparisons
# =========================================================

resultsNames(dds2)
# =========================================================
# GSEA Analysis: Healthy vs RA in M1-like LPS
# =========================================================

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)

# =========================================================
# Extract DESeq2 Results
# =========================================================
# RA vs Healthy within M1-like LPS

res_gsea <- results(
  dds2,
  contrast = c(
    "Disease_Group",
    "rheumatoid arthritis_M1-like_LPS",
    "healthy_M1-like_LPS"
  )
)

# Convert to dataframe
res_gsea <- as.data.frame(res_gsea)

# Remove NA values
res_gsea <- na.omit(res_gsea)

# Add gene IDs
res_gsea$ENSEMBL <- rownames(res_gsea)

# Remove ENSEMBL Version Numbers

res_gsea$ENSEMBL <- gsub(
  "\\..*",
  "",
  res_gsea$ENSEMBL
)

# Convert ENSEMBL → ENTREZ

gene_map <- bitr(
  res_gsea$ENSEMBL,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Merge mapping
res_gsea <- merge(
  res_gsea,
  gene_map,
  by.x = "ENSEMBL",
  by.y = "ENSEMBL"
)

# Create Ranked Gene List
# Ranking by log2FoldChange

gene_list <- res_gsea$log2FoldChange

names(gene_list) <- res_gsea$ENTREZID

# Sort decreasing
gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

# Run GSEA GO Analysis
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View Results

head(as.data.frame(gsea_go))

# Save Results
write.csv(
  as.data.frame(gsea_go),
  "GSEA_GO_RA_vs_Healthy_M1likeLPS.csv",
  row.names = FALSE
)



# KEGG GSEA

gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Save KEGG Results

write.csv(
  as.data.frame(gsea_kegg),
  "GSEA_KEGG_RA_vs_Healthy_M1likeLPS.csv",
  row.names = FALSE
)

# KEGG Visualization

# GO GSEA dotplot
dotplot(
  gsea_go,
  showCategory = 15,
  title = "GSEA GO: RA vs Healthy in M1-like LPS"
)

# Enrichment curve
gseaplot2(
  gsea_go,
  geneSetID = 1,
  title = gsea_go$Description[1]
)

# Ridgeplot
ridgeplot(
  gsea_go,
  showCategory = 15
))



# =========================
# Finish
# =========================

