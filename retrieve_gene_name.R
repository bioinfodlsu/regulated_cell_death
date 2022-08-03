data = read.csv("output/DE_Results_ordered_gene_names.csv")
data
colnames(data)[1] <- "ensembl_gene_id"
head(data)
library('biomaRt')

hsapiens = useMart("ensembl",
                   dataset="hsapiens_gene_ensembl")

hsapiens_infos <- getBM(attributes=c('ensembl_gene_id',
                                     'external_gene_name'),
                        mart = hsapiens)
hsapiens_infos
merge_infos <- merge(x = data,
                     y = hsapiens_infos,
                     by = "ensembl_gene_id",
                     all.x = TRUE)

head(merge_infos)
merge_infos<- merge_infos[order(merge_infos$padj),]
head(merge_infos)
write.csv(merge_infos, file = "output/gene_names.csv")
