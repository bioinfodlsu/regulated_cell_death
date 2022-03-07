#gene set enrichment analysis

loadLibraries <- function()
{
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggridges)
  library(org.Hs.eg.db)
  require(DOSE)
  
}

loadDeseqFile <- function(file)
{
  df = read.csv(file, header=TRUE)
  gene_list <- df$log2FoldChange
  names(gene_list) <- df$X
  retList <- list("gene_list"= gene_list, "df"=df)
  return(retList)
}

performGOAnalysis <-function(gene_list, minGSSize= 3, maxGSSize = 800, pvalueCutoff = 0.05)
{
  

  go_gene_list<-na.omit(gene_list)
  go_gene_list = sort(go_gene_list, decreasing = TRUE)
  
  set.seed(42)
  
  gse <- gseGO(geneList=go_gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               minGSSize = minGSSize, 
               maxGSSize = maxGSSize, 
               pvalueCutoff = pvalueCutoff,
               seed = TRUE,
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH")
  return (gse)
}

performKEGGAnalysis <-function(gene_list, df, minGSSize= 3, maxGSSize = 800, pvalueCutoff = 0.05)
{
  organism = "org.Hs.eg.db"
  ids<-bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  df2 = df[df$X %in% dedup_ids$ENSEMBL,]
  df2$Y = dedup_ids$ENTREZID
  kegg_gene_list <- df2$log2FoldChange
  names(kegg_gene_list) <- df2$Y
  kegg_gene_list<-na.omit(kegg_gene_list)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  kegg_organism = "hsa"
  set.seed(42)
  
  kegg <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 minGSSize    = minGSSize,
                 maxGSSize    = maxGSSize,
                 pvalueCutoff = pvalueCutoff,
                 seed = TRUE,
                 pAdjustMethod = "BH",
                 keyType       = "ncbi-geneid")
  return(kegg)
}

createdotPlot <- function(data, numTerms, title)
{
  if (missing(numTerms))
  {
    dotplot(data,title = title, split=".sign") + facet_grid(.~.sign)
  }
  else
  {
    dotplot(data, showCategory=numTerms, title = title, split=".sign") + facet_grid(.~.sign)
  }
}

createRidgePlot <- function(data, xlabel)
{
  ridgeplot(data) + labs(x = xlabel)
}

createEM <- function(data, numTerms)
{
  x <- pairwise_termsim(data) 
  if (missing(numTerms))
  {
    emapplot(x)
  }
  else
  {
    emapplot(x, showCategory=numTerms)
  }
}

createPathView <- function (pathwayID)
{
  library(pathview)
  dme <- pathview(gene.data=kegg_gene_list, pathway.id=pathwayID, species = kegg_organism)
  dme <- pathview(gene.data=kegg_gene_list, pathway.id=pathwayID, species = kegg_organism, kegg.native = F)
  
}

  #sample code for GSEA
  loadLibraries()
  data = loadDeseqFile("output/DE_Results_ordered_gene_names.csv")
  go = performGOAnalysis(gene_list=data$gene_list, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05)

  createdotPlot(go, 10, NULL)
  createEM(go,10)
  createRidgePlot(go, "enrichment distribution")
  
  kegg = performKEGGAnalysis(gene_list = data$gene_list,
            df = data$df,
            minGSSize    = 3,
            maxGSSize    = 800,
            pvalueCutoff = 0.05)
  
  createdotPlot(kegg, 10, "KEGG Enriched Pathways")
  createEM(kegg,10)
  createRidgePlot(kegg, "enrichment distribution")

  #pathview
  createPathView("hsa00982")
