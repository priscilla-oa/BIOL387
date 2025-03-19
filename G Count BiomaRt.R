
install.packages("openxlsx")
library(biomaRt)
library(stringr)
library(readxl)
library(openxlsx)  


human <- useMart(biomart = "ensembl", 
                 dataset = "hsapiens_gene_ensembl", 
                 verbose = TRUE, 
                 host = "https://dec2021.archive.ensembl.org")
data <- read_excel("/Users/priscilla/Documents/Dry Lab/Gene_id.xlsx") 
head(data)
results <- list()
for (gene in data$gene_name) { 
  gene_info <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                     filters = "hgnc_symbol",
                     values = gene,
                     mart = human)
  if (nrow(gene_info) == 0) {
    next  
  }
  transcript <- getSequence(id = gene_info$ensembl_gene_id, 
                            type = 'ensembl_gene_id', 
                            seqType = "gene_exon_intron", 
                            mart = human)
  if (is.null(transcript) || length(transcript) == 0) {
    next
  }
  transcript <- as.character(transcript)
  g_count <- sum(str_count(transcript, "G"))
  results[[gene]] <- data.frame(gene_name = gene, g_count = g_count)
}
final_results <- do.call(rbind, results)
print(final_results)
write.xlsx(final_results, "/Users/priscilla/Documents/Dry Lab/Ensembl_Gene_G_Counts.xlsx")
