library(tidyverse)

df <- read_delim("~/qbb2024-answers/d4-morning/genes_tissues_expression.tsv")

view(df)



df_tissue_geneid <- df %>%
  dplyr::mutate(Tissue_Gene=paste0(Tissue, " ", Gene_ID))

view(df_tissue_geneid)


df_log2 <- df_tissue_geneid %>%
  dplyr::mutate(log2_expression = log2(Expression + 1.0))

view(df_log2)


ggplot(data = df_log2, mapping = aes(x = Tissue_Gene, y = log2_expression)) + 
  geom_violin() +
  coord_flip() +
  labs(
    title = "Expression variability of key genes in various tissues",
    x = "Tissue and Gene of Interest",
    y = "Gene Expression (log2)"
  )


# Given the high expression and high tissue specificity observed for these 
    # genes, I would have expected these genes of interest to be essential 
    # genes. As such, I was expecting each gene to have incredibly low 
    # expression variability within a population, ensuring that all members 
    # of the population have relatively equal levels of the gene products. 
    # Surprisingly, the violin plots don't support that hypothesis, and 
    # instead show that most of the genes of interest have remarkably high 
    # expression variability.

# Despite the observations described above, expression of highly 
    # tissue-specific genes in the pancreas is consistent within the 
    # population (very low variability). This suggests that those pancreatic 
    # genes are essential for survival, meaning everyone in a given 
    # population should be expected to have relatively equal levels of the 
    # gene products. 

