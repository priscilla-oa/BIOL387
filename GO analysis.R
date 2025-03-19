

#Code for GO Molecular Function and Reactome
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
gene_list <- c("MYH7", "MYBPC3", "FLNC", "DSP", "TTN", "TCAP", "LZTR1", "TNNI3", "COL3A1", 
               "APOB", "DOLK", "SCO2", "FBN1", "EMD", "LDLR")
gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_result <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05)
reactome_result <- enrichPathway(gene = gene_ids$ENTREZID, organism = "human", pvalueCutoff = 0.05)
head(go_result)
head(reactome_result)
barplot(go_result, showCategory = 10)
barplot(reactome_result, showCategory = 10)

#to Visualise genes on Reactome Plot
top_15_genes <- c("MYH7", "MYBPC3", "FLNC", "DSP", "TTN", "TCAP", "LZTR1", 
                  "TNNI3", "COL3A1", "APOB", "DOLK", "SCO2", "FBN1", "EMD", "LDLR")

reactome_df <- as.data.frame(reactome_result)
reactome_df$GeneListSymbols <- sapply(reactome_df$geneID, function(gene_id) {
  genes <- unlist(strsplit(gene_id, "/"))  
  symbols <- bitr(genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db) 
  if (length(symbols$SYMBOL) > 0) {
    paste(symbols$SYMBOL, collapse = ", ")  
  } else {
    NA 
  }
})
head(reactome_df[, c("Description", "GeneListSymbols")])
library(ggplot2)
library(stringr)
reactome_df$DescriptionWrapped <- str_wrap(reactome_df$Description, width = 40)  
ggplot(reactome_df, aes(x = reorder(DescriptionWrapped, -pvalue), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = GeneListSymbols), hjust = 1.1, color = "black", size = 3) +
  labs(x = "Pathway", y = "-log10(P-value)", title = "Top Reactome Pathways Enriched with Genes") +
  coord_flip() +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 10, hjust = 0.5),  
    plot.margin = margin(10, 10, 20, 10) 
  )
