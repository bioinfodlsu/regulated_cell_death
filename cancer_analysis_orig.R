library("vsn")
directory <- "/Users/jenny/Desktop/Cancer Project/Cancer Project/data"
outdir <- "Users/jenny/Desktop/Cancer Project/Cancer Project/output"
directory
list <-list.files(directory)
list
sampleFiles <- grep("normal",list.files(directory),value=TRUE)
head(sampleFiles)
sampleCondition <- sub("(.*normal).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

#filter to keep only rows that have 10 reads total. before: 60483 x 519 after: 53996 x 519
keep <- rowSums(counts(ddsHTSeq)) >= 10  
ddsHTSeq <- ddsHTSeq[keep,]
head(ddsHTSeq,10)

#for visualization or clustering – it might be useful to work with transformed versions of the count data
#if you have many samples (e.g. 100s), the rlog function might take too long, and so the vst function will be a faster choice
vsd <- vst(ddsHTSeq, blind=FALSE)

#plot the standard deviation of the transformed data, across samples
#against the mean, using  the variance stabilizing transformation
#the vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions
sdplot <-meanSdPlot(assay(vsd))
library("ggplot2")
sdplot$gg <- sdplot$gg +
  ggtitle(label="Mean SD of transformed read counts") +
  ylab("standard deviation")
print(sdplot$gg)

#saving plot not allowed on rstudio
#savePlot(filename = "Mean SD", type = "png",)

#data quality assessment
#heatmap of count matrix
#library("pheatmap")
#library(grid)
#select <- order(rowMeans(counts(ddsHTSeq,normalized=TRUE)),
#                decreasing=TRUE)[1:20]
#rownames(df) <- colnames(ddsHTSeq)
#heatmap<-pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
#print(heatmap)
#str(df)
#str(select)

#PCA Plot
pcaData<-plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#to further check outliers
boxplot(log10(assays(ddsHTSeq$condition=="normal")[["cooks"]]), range=0, las=2)

#level represents the control group or which level you want to compare against
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = c("normal","notnormal")) 
#differential expression testing
ddsHTSeq <- DESeq(ddsHTSeq)
res <- results(ddsHTSeq)
res
write.csv(res, file = "output/DE_Results.csv")

#Log fold change shrinkage for visualization and ranking
resultsNames(ddsHTSeq)
resLFC <- lfcShrink(ddsHTSeq, coef="condition_notnormal_vs_normal", type="apeglm")
write.csv(resLFC, file = "output/DE_Results_lfc.cs")

#order our results table by the smallest p value:
resOrdered <- resLFC[order(resLFC$padj),]
resOrdered
write.csv(resOrdered, file = "output/DE_Results_ordered.csv")

#By default the argument alpha is set to 0.1
res05 <- results(ddsHTSeq, alpha=0.05)
sink("output/log.txt")
print(summary(res05))
sink()
#read counts per gene
#You can select the gene to plot by rowname or by numeric index
#this sample gets the gene which had the smallest p value 
plotCounts(ddsHTSeq, gene=which.min(resOrdered$padj), intgroup="condition")

resOrdered
head(resOrdered,10)
#get top 10 genes
colnames(resOrdered)
goi <- rownames(resOrdered)[1:10]
goi
for(i in 1:10) {
  d <- plotCounts(ddsHTSeq, gene=goi[i], intgroup="condition", returnData = )
}
plotCounts(ddsHTSeq, gene=which.min(res$padj), intgroup="condition")