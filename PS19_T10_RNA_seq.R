library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggview)
library(writexl)
ggview <- ggview:::ggview


# Load data
countData <- as.matrix(read.csv("T10.csv", row.names = "Gene"))

# Remove genes with low expression levels
countData <- countData[rowMeans(countData) > 10,]

condition <- factor(c(rep("WT",6), rep("WT_T10_kd",6), rep("PS19",5), rep("PS19_T10_kd",4)))
colData <- data.frame(row.names = colnames(countData), condition)

# Differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~ condition)
dds1 <- DESeq(dds, fitType = "parametric", minReplicatesForReplace = 7, parallel = FALSE)

# Extract results using the results() function
res <- results(dds1, contrast = c("condition", "PS19_T10_kd", "PS19"))

# Convert the format of 'res' to a data frame using data.frame()
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

# Sort sequentially by padj value and then by log2FoldChange value
DEG <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

# Remove missing values
DEG1 <- na.omit(DEG)


# Data normalization (variance-stabilizing transformation)
vsd <- vst(dds1, blind = FALSE)

# Generate data for PCA plotting
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the variance explained by each principal component
percentVar <- round(100 * attr(pcaData, "percentVar"))

P <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 8) +
  theme_bw() +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_vline(xintercept=0, linetype = "dashed")+
  theme(legend.position = "none",
        axis.title = element_text(vjust = 0, size = 25),
        axis.text = element_text(color = "black", size = 25),
        axis.ticks = element_line(color = "black", linewidth = 2),
        panel.border = element_rect(colour = "black", linewidth = 3))
ggview(P, width = 10, height = 9)
ggsave(plot = P, filename = "PCA-T10.pdf", width = 10, height = 9)


# Identify differentially expressed genes
DEG1_up <- DEG1[which(DEG1$log2FoldChange >= 0.4 & DEG1$padj < 0.05),]
DEG1_down <- DEG1[which(DEG1$log2FoldChange <= -0.4 & DEG1$padj < 0.05),]
DEG1_total <- rbind(DEG1_up[1:50, ], DEG1_down[1:70, ])

# Locate differentially expressed genes in the original expression matrix
df <- countData[intersect(rownames(countData), rownames(DEG1_total)),]

# Convert data frame to matrix format
df1 <- as.matrix(df)

# Plot heatmap
P1 <- pheatmap(df1,
               scale = "row",                                           # Row-wise normalization
               border="white",                                          # Set border color to white
               cluster_cols = F,                                        # Disable column clustering
               show_rownames = F,                                       # Hide row names
               show_colnames = F,                                       # Hide column names
               fontsize = 24,                                           # Set font size
               color = colorRampPalette(c("blue","white","red"))(100),  # Custom color palette
) 
ggview(P1, width = 10, height = 9)
ggsave(plot = P1, filename = "Heatmap-T10.pdf", width = 10, height = 9)


DEG1_up_with_rownames <- cbind(gene = rownames(DEG1_up), DEG1_up)
write_xlsx(DEG1_up_with_rownames, "DEG_UP.xlsx")

# KEGG enrichment analysis of differentially expressed genes using DAVID
# Plot KEGG image
UP_KEGG <- read.csv(file = "KEGG_UP.csv")
UP_KEGGBP <- subset(UP_KEGG, subset = (Ontology == "KEGG_PATHWAY"))[1:11,]
UP_KEGGBP$Description <- factor(UP_KEGGBP$Description, levels = rev(UP_KEGGBP$Description))

P2 <-ggplot(data = UP_KEGGBP,
            aes(x = Ratio, y = reorder(Description, Count))) +
  scale_size_continuous(range = c(2, 10)) + 
  geom_point(aes(size = Count, color = log10padj)) +
  theme_bw() +
  scale_colour_gradient(low = "yellow",high = "red") +
  scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
  labs(x = "GeneRatio",y = "",
       color = expression(-log10padj),size = "Count") +
  guides(colour = guide_colorbar(barwidth = 1, barheight = 9, order = 1),
         size = guide_legend(keywidth = 1, keyheight = 1.7, order = 2)) +
  theme(axis.title = element_text(color = "black", size = 24),
        axis.text = element_text(color = "black", size = 22),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 25,hjust = 0.5),
        plot.margin = margin(l = -0.7, unit = "cm"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank())
ggview(P2, width = 11, height = 9)
ggsave(plot = P2, filename = "KEGG-T10.pdf", width = 11, height = 9)

dev.off()
