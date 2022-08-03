#gene set enrichment analysis
#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# reading in data from deseq2
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(ggridges)
library(org.Hs.eg.db)
organism = "org.Hs.eg.db"

df = read.csv("output/DE_Results_ordered_gene_names.csv", header=TRUE)
df
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X


#KEGG Enrichment analysis
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted #40% not mapped if entrezid 65% if uniprot
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% dedup_ids$ENSEMBL,]
df2
notmapped = df[!(df$X %in% dedup_ids$ENSEMBL),]
notmapped
write.csv(notmapped, file = "output/ensembl_nomapping.csv")
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$SYMBOL

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)



kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_gene_list
y <- GSEA(kegg_gene_list, TERM2GENE=gmt)
y
write.csv(y, file = "output/ferroptosis_gsea.csv")
kegg_gene_list
gmt <- read.gmt("WP_FERROPTOSIS.v7.5.1.gmt")
wp <- gseWP(kegg_gene_list, organism = "Homo sapiens",TERM2GENE=gmt)
wp_ordered <- wp[order(wp$p.adjust),]
write.csv(wp, file = "output/wikipathways_result_ordered.csv")


dotplot(wp, showCategory = 10,label_format=100, title = "WikiPathways Enriched Pathways" , split=".sign") + facet_grid(.~.sign)


x2 <- pairwise_termsim(wp) 
emapplot(x2,label_format=100,cex_label_category=1,layout = "nicely",showCategory = 20, node_scale = NULL,
         line_scale = NULL,
         min_edge = 0.2)
ridgeplot(wp,label_format=100,showCategory = 20)



gmt
y <- GSEA(kegg_gene_list, TERM2GENE=gmt,pvalueCutoff = 1)
