# Bulk-seq
For KEGG pathway analysis from DEseq2 results of control vs. experimental
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # For human genes
BiocManager::install("pathview")      # For pathway visualization (optional)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Read the CSV file
gene_data <- read.csv("C://BulkRNA//res.csv")

# Ensure correct column names
head(gene_data)


# Filter upregulated genes
upregulated_genes <- gene_data %>% 
  filter(log2FoldChange > 1) %>% 
  pull(X) # X is ENTREZID column head


# Filter downregulated genes
downregulated_genes <- gene_data %>% 
  filter(log2FoldChange < -1) %>% 
  pull(X)

# Check the number of genes in each category
cat("Number of upregulated genes:", length(upregulated_genes), "\n")
cat("Number of downregulated genes:", length(downregulated_genes), "\n")


# KEGG analysis for upregulated genes
kegg_up <- enrichKEGG(
  gene = upregulated_genes,
  organism = 'hsa',    # 'hsa' is for Homo sapiens
  keyType = 'kegg',    # Key type is KEGG IDs
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2
)

# KEGG analysis for downregulated genes
kegg_down <- enrichKEGG(
  gene = downregulated_genes,
  organism = 'hsa',
  keyType = 'kegg',
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)


# Bar plot for upregulated pathways
barplot(kegg_up, showCategory = 10, title = "Upregulated KEGG Pathways for HSA data") 

# Bar plot for downregulated pathways
barplot(kegg_down, showCategory = 10, title = "Downregulated KEGG Pathways for HSA data")

# Dot plot for better visualization
dotplot(kegg_up, showCategory = 10, title = "Upregulated KEGG Pathways")
dotplot(kegg_down, showCategory = 10, title = "Downregulated KEGG Pathways")


# Save upregulated results
write.csv(as.data.frame(kegg_up), "C://BulkRNA//KEGGres//kegg_upregulated.csv", row.names = FALSE)

# Save downregulated results
write.csv(as.data.frame(kegg_down), "C://BulkRNA//KEGGres//kegg_downregulated.csv", row.names = FALSE)


library(pathview)

# Visualize the top enriched pathway for upregulated genes
pathview(gene.data = upregulated_genes, pathway.id = kegg_up@result$ID[1:10], species = "hsa")

# Visualize the top enriched pathway for downregulated genes
pathview(gene.data = downregulated_genes, pathway.id = kegg_down@result$ID[1:10], species = "hsa")
