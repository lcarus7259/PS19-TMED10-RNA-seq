# Differential expression analysis of RNA-seq data using DESeq2

This repository contains the complete R code to reproduce the differential expression analysis, PCA visualization, heatmap, and KEGG enrichment plots described in the manuscript. The analysis was performed using **DESeq2** (Love et al., 2014) on RNA-seq count data from four experimental groups.

## Requirements

- **R** version 4.3 or higher
- Required R packages (install if missing):
  ```r
  install.packages(c("ggplot2", "pheatmap", "ggrepel", "dplyr", "stringr", "writexl"))
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("DESeq2", "ggview"))
  ```

## Input data

Place the following file in the same directory as the R script:

- **`T10.csv`**  
  A raw count matrix (genes as rows, samples as columns).  
  - First column must be named `"Gene"` and contains gene identifiers (e.g., Ensembl IDs or gene symbols).  
  - All other columns contain integer read counts for each sample.  

The sample condition (grouping) is defined inside the script as:
```r
condition <- factor(c(rep("WT",6), rep("WT_T10_kd",6), rep("PS19",5), rep("PS19_T10_kd",4)))
```
Make sure the order of columns in `T10.csv` matches this condition vector (i.e., first 6 samples WT, next 6 WT_T10_kd, next 5 PS19, last 4 PS19_T10_kd).

## Running the analysis

Open R or RStudio and run the entire script:
   ```r
   source("PS19_T10_RNA_seq.R")
   ```

The script performs the following steps automatically:

- Loads and filters low‑expression genes (row mean > 10)
- Runs DESeq2 differential expression analysis (comparison: `PS19_T10_kd` vs `PS19`)
- Orders results by adjusted p‑value and log2 fold change
- Generates a PCA plot (`PCA-T10.pdf`)
- Identifies differentially expressed genes (|log2FC| ≥ 0.4 and padj < 0.05)
- Produces a heatmap using the top 50 up‑regulated and top 70 down‑regulated genes (`Heatmap-T10.pdf`)
- Exports up‑regulated genes to `DEG_UP.xlsx` (used for KEGG enrichment)

### KEGG enrichment plotting (additional step)

The script reads a pre‑computed KEGG enrichment file named **`KEGG_UP.csv`** that you must generate separately using **DAVID**.  

Expected columns: `Ontology`, `Description`, `Ratio`, `Count`, `log10padj`.  

To obtain this file:

- Take the gene list from `DEG_UP.xlsx` (gene column)
- Submit it to DAVID (https://davidbioinformatics.nih.gov/) for KEGG pathway analysis, **all genes detected in experiment** (i.e., the gene list from `T10.csv`) is used as the background
- Download the result table and save as `KEGG_UP.csv` in the same folder
- Ensure the column `Ontology` contains "KEGG_PATHWAY" for the rows you wish to plot
- `Ratio` was calculated as `Count` divided by the total number of genes
- The column `log10padj` actually represents **`-log10(padj)`** (the negative base‑10 logarithm of the adjusted p‑value), where larger values indicate more significant enrichment

After placing `KEGG_UP.csv`, re‑run the script:
```r
source("PS19_T10_RNA_seq.R")
```

The script then produces `KEGG-T10.pdf` (bubble plot of top 11 KEGG pathways).

## Output files

| File | Description |
|------|-------------|
| `PCA-T10.pdf` | Principal component analysis (PC1 vs PC2) colored by condition. |
| `Heatmap-T10.pdf` | Heatmap of top differentially expressed genes (row‑scaled). |
| `KEGG-T10.pdf` | Bubble plot for KEGG enrichment (GeneRatio vs Count, coloured by significance). |
| `DEG_UP.xlsx` | Excel file with up‑regulated genes (log2FC ≥ 0.4, padj < 0.05). |

## Reproducibility

To ensure identical results, run the following to see your R environment:

```r
sessionInfo()
```

Example output from the original analysis:
```
R version 4.3.0 (2023-04-21 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26200)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base    

other attached packages:
 [1] writexl_1.5.4               ggview_0.2.2                stringr_1.5.0              
 [4] dplyr_1.1.2                 ggrepel_0.9.3               pheatmap_1.0.12            
 [7] ggplot2_3.4.3               DESeq2_1.40.2               SummarizedExperiment_1.30.2
[10] Biobase_2.60.0              MatrixGenerics_1.12.3       matrixStats_1.0.0          
[13] GenomicRanges_1.52.0        GenomeInfoDb_1.36.1         IRanges_2.34.1             
[16] S4Vectors_0.38.1            BiocGenerics_0.46.0
```

**Note:** The exact output may vary slightly with newer R/package versions, but the overall biological conclusions should remain unchanged.

## License

Copyright © 2026 Authors.

This repository is provided solely for academic reproducibility.
No commercial use, redistribution, or modification is permitted without prior permission from the authors.

## Citation

If you use this analysis in your work, please cite:
- Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12):550.
