### Description

DESeq2<sup>[1]</sup> is an open source [R package](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) for differential analysis of count data, using shrinkage estimation for dispersions and fold changes to improve stability and interpretability of estimates. This enables a more quantitative analysis focused on the strength rather than the mere presence of differential expression.

DESeq2 offers a comprehensive and general solution for gene-level analysis of RNA-seq data. Shrinkage estimators substantially improve the stability and reproducibility of analysis results compared to maximum-likelihood-based solutions. Empirical Bayes priors provide automatic control of the amount of shrinkage based on the amount of information for the estimated quantity available in the data. This allows DESeq2 to offer consistent performance over a large range of data types and makes it applicable for small studies with few replicates as well as for large observational studies.

### Installation  

First, we need to install `BiocManager`:  

    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
   
Then we just call  

    BiocManager::install("DESeq2")
    library(DESeq2)

### References

[1] Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550 (2014). DOI:10.1186/s13059-014-0550-8
