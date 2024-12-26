#Read count data
library(readr)
counts1.a <- read.table('/Users/sruthi/Downloads/FILE1-COUNTS1.txt', sep="\t", header=FALSE, stringsAsFactors=FALSE)
head(counts1.a)
counts1 <- as.matrix(counts1.a[-1,-1])
mode(counts1) <- "integer"   # Convert to integer
rownames(counts1) <- counts1.a[-1,1]     # Set rownames to geneIDs
colnames(counts1) <- counts1.a[1,-1]
ncol(counts1)
#Sample information
colData1.a <- read.csv('/Users/sruthi/Downloads/FILE2-COLDATA1.csv')
nrow(colData1.a)
colData1 <- as.matrix(colData1.a[-1,-1])
rownames(colData1) <- colData1.a[-1,1]
colnames(colData1) <- colData1.a[1,-1]
View(colData1)
# check if the data is ready for DESeq2 analysis
library(DESeq2)
all(colnames(counts1) %in% rownames(colData1))
all(colnames(counts1) == rownames(colData1)) 
#DESeq2 Analysis
dds1.0 <- DESeqDataSetFromMatrix(countData = counts1, 
                                 colData = colData1, 
                                 design = ~ condition)
keep <- rowSums(counts(dds1.0)) >= 10
dds1 <- dds1.0[keep,]
dds1
# set the factor level
dds1$condition <- relevel(dds1.0$condition, ref = " Control")
levels(dds1$condition)
#Step 3: Run DESeq
dds1 <- DESeq(dds1)
res <- results(dds1)
summary(res)
normalizedcounts1 <- counts(dds1, normalized = TRUE)
head(normalizedcounts1)
write.csv(normalizedcounts1, 'normalizedcounts1.csv')
#Alpha = 0.05
res1 <- results(dds1, alpha = 0.05)
res1
summary(res1)
resOrdered1 <- res1[order(res1$padj),]
resOrdered1
#plotting the genes with increasing the size of the dots on the graph by cex=0.7
plotMA(res1,cex =0.3, ylim=c(-10,10))
abline (h=c(-1,1), col="red", lwd=3)
#Dispersion plot
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
library(apeglm)
resultsNames(dds1)
res1LFC <- lfcShrink(dds1, coef = "condition_.sample_vs_.Control", type='apeglm')
plotMA(res1LFC, cex = 0.5, ylim=c(-10,10))
plotDispEsts(dds1, main = "Dispersionplot1")

#PCA analysis
install.packages("ggplot2")
library(ggplot2)
library(DESeq2)
rld1 <- rlog(dds1, blind = FALSE)
head(assay(rld1))
hist(assay(rld1))
PCAA1 <- plotPCA(rld1,intgroup='condition')
PCAA1 + geom_text(aes(label = name), size = 2)+ggtitle('PCA plot1')

library(ggplot2)
padj.cutoff <- 0.05
lfc.cutoff <- 2
threshold.1 <- res1$padj < padj.cutoff & abs(res1$log2FoldChange) > lfc.cutoff
res1$threshold_S1 <- threshold.1
sig.S1 <- data.frame(subset(res1, threshold_S1==TRUE))
sig.S1_ordered <- sig.S1[order(sig.S1$padj), ]
top20_sig.S1_genes <- rownames(sig.S1_ordered)[1:20]
View(top20_sig.S1_genes)
normalizedcounts1 <- counts(dds1, normalized = TRUE)
avg_normalizedcounts_sample1 <- rowMeans(normalizedcounts1)
top20_sigs1_norm <- as.data.frame(normalizedcounts1[top20_sig.S1_genes, ])
View(top20_sigs1_norm)
write.csv(top20_sigs1_norm, 'top20_sigs1_norm')

######################################data2
#Read the data
library(readr)
Counts_data2 <- read.table('~/Downloads/FILE3-COUNTS2.csv')
counts2 <- as.matrix(Counts_data2[-1,-1])
mode(counts2) <- "integer"   # Convert to integer
rownames(counts2) <- Counts_data2[-1,1]     # Set rownames to geneIDs
colnames(counts2) <- Counts_data2[1,-1]
ncol(counts2)
View(counts2)
#Sample information
colData2.0 <- read.csv('~/Downloads/FILE4-COLDATA2 (1).csv')
View(colData2.0)
colData2 <- as.matrix(colData2.0[-1,-1])
rownames(colData2) <- colData2.0[-1,1]
colnames(colData2) <- colData2.0[1,-1]
View(colData2)
#Check if the data2 is ready for DESeq2
all(colnames(counts2) %in% rownames(colData2))
all(colnames(counts2) == rownames(colData2)) 
library(DESeq2)
#DESeq2 Analysis
dds2.0 <- DESeqDataSetFromMatrix(countData = counts2, 
                                 colData = colData2, 
                                 design = ~ Condition)
keep <- rowSums(counts(dds2.0)) >= 10
dds2 <- dds2.0[keep,]
dds2
# set the factor level
dds2$Condition <- relevel(dds2$ Condition, ref = "Control")
levels(dds2$Condition)
dds2 <- DESeq(dds2)
res2<- results(dds2, alpha = 0.05)
res2
summary(res2)
plotMA(res2)
normalizedcounts2 <- counts(dds2, normalized = TRUE)
head(normalizedcounts2)
write.csv(normalizedcounts2, 'normalizedcounts2.csv')
resOrdered2 <- res2[order(res2$padj),]
resOrdered2
#plotting the genes with increasing the size of the dots on the graph by cex=0.7
plotMA(res2,cex =0.3, ylim=c(-10,10))
abline (h=c(-1,1), col="red", lwd=3)
#Dispersion plot
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
library(apeglm)
resultsNames(dds2)
res2LFC <- lfcShrink(dds2, coef = "Condition_Treated_vs_Control", type='apeglm')
plotMA(res2LFC, cex = 0.5, ylim=c(-10,10))
plotDispEsts(dds2, main = "Dispersionplot2")
#PCA 
suppressPackageStartupMessages(library(ggplot2))
install.packages("ggplot2")
library(ggplot2)
library(DESeq2)
rld2 <- rlog(dds2, blind = FALSE)
head(assay(rld2))
hist(assay(rld2))
PCAA2 <- plotPCA(rld2,intgroup='Condition')
PCAA2 + geom_text(aes(label = name), size = 2)+ggtitle('PCA plot2')
#Acquiring top 20 highly expressed genes
padj.cutoff <- 0.05
lfc.cutoff <- 2
threshold.2 <- res2$padj < padj.cutoff & abs(res2$log2FoldChange) > lfc.cutoff
res2$threshold_S2 <- threshold.2
sig.S2 <- data.frame(subset(res2, threshold_S2==TRUE))
sig.S2_ordered <- sig.S2[order(sig.S2$padj), ]
top20_sig.S2_genes <- rownames(sig.S2_ordered)[1:20]
View(top20_sig.S2_genes)
normalizedcounts2 <- counts(dds2, normalized = TRUE)
avg_normalizedcounts_sample2 <- rowMeans(normalizedcounts2)
top20_sigs2_norm <- as.data.frame(normalizedcounts2[top20_sig.S2_genes, ])
View(top20_sigs2_norm)
write.csv(top20_sigs2_norm, 'top20_sigs2_norm')
View(top20_sigs2_norm)


