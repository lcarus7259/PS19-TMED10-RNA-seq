library(ggplot2)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(stringr)

# 读入数据，注意设置工作路径
countData <- as.matrix(read.csv("T10.csv",row.names="Gene"))
# 去除表达量过低的基因
countData <- countData[rowMeans(countData)>10,]

condition <- factor(c(rep("WT",6),rep("WT_T10_kd",6),rep("PS19",5),rep("PS19_T10_kd",4)))
colData <- data.frame(row.names=colnames(countData), condition)
dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~ condition)

# 差异表达分析
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
# 将结果用result()函数来获取
res <- results(dds1, contrast=c("condition","PS19","WT"))
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照padj值log2FoldChange值进行排序
DEG <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
# 去除缺失值
DEG_deseq2 <- na.omit(DEG)


#绘制KEGG图像
UP_KEGG <- read.csv(file = "C:/Users/Icarus/Desktop/KEGG_UP_SHORT.csv")
UP_KEGGBP <- subset(UP_KEGG, subset = (Ontology == "KEGG"))[1:8,]
UP_KEGGBP$Description <- factor(UP_KEGGBP$Description, levels = rev(UP_KEGGBP$Description))
pdf(file = 'KEGG-T10.pdf', width = 11, height = 9)
ggplot(data = UP_KEGGBP,
  aes(x = Ratio, y = reorder(Description, Count))) +
  scale_size_continuous(range = c(2, 10)) + 
  geom_point(aes(size = Count,color = log10padj)) +
  theme_bw() +
  scale_colour_gradient(low = "blue",high = "red") +
  scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
  labs(x = "GeneRatio",y = "",title = "PS19; T10fl/+ vs PS19",
       color = expression(-log10padj),size = "Count") +
  guides(colour = guide_colorbar(barwidth = 1, barheight = 9),
         size = guide_legend(keywidth = 1, keyheight = 1.7)) +
  theme(axis.title = element_text(color = "black", size = 25),
        axis.text = element_text(color = "black", size = 25),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 25,hjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid.minor = element_blank())
dev.off()  

# 数据标准化
vsd <- vst(dds1, blind=FALSE)
#生成pca作图数据
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
#计算主成分的方差贡献率
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = 'PCA-T10.pdf', width = 10, height = 9)
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 8) +
  theme_bw() +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_vline(xintercept=0, linetype = "dashed")+
  theme(legend.position = "none",
        axis.title = element_text(vjust = 0, size = 25),
        axis.text = element_text(color = "black", size = 25),
        axis.ticks = element_line(color = "black", size = 2),
        panel.border = element_rect(colour = "black", linewidth = 3))
dev.off()


