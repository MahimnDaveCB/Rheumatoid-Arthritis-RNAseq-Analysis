# Transcriptomic Analysis of peripheral Monocytes Under Distinct Activation Conditions isolated from RA patients

## Overview

This project performs differential gene expression and transcriptomic profiling of peripheral blood CD14+ monocytes isolated from Rheumatoid Arthritis (RA) patients and matched healthy donors under multiple ex vivo activation conditions using RNA-seq data.

The analysis investigates disease-associated transcriptional alterations across macrophage differentiation states and inflammatory stimulation conditions using DESeq2-based RNA-seq workflows.

# Dataset Information
- **Title:** Transcriptomic comparison of RA vs healthy monocytes undergoing different ex vivo activation conditions
- **Organism:** Homo sapiens
- **Experiment Type:** Expression profiling by high throughput sequencing
- **Platform:** RNA-seq
- **Data Source:** GEO Dataset (GSE328878)
- **Data Contributors:** Teoh S, Börsch A, Müller-Durovic B

## Repository Structure

```text
Rheumatoid-Arthritis-RNAseq-Analysis/
│
├── data/
│   ├── raw_counts.csv
│   └── design.csv
│
├── figures/
│   ├── DEGs Heatmap.png
│   ├── Dispersion estimates.png
│   ├── Frequencies of padj values.png
│   ├── GO Biological Process.png
│   ├── KEGG pathway enrichment.png
│   ├── PCA Plot.png
│   ├── PPI interaction network.png
│   └── Volcano Plot M2-like LPS vs Vehicle.png
│
├── results/
│   ├── deseq.result.all.tsv
│   ├── deseq_deg.tsv
│   ├── normalized_counts.csv
│   ├── GO_enrichment.tsv
│   └── KEGG_enrichment.tsv
│
├── scripts/
│   └── RNAseq_DESeq2_pipeline.R
│
└── README.md
```
