library(tidyverse)
library(ggthemes)

# read in the snp_counts file generated in exercise 2.1 and save as a new dataframe
snp_enrichment <- read_tsv("~/qbb2024-answers/week2/snp_counts.txt")

# add a column to the snp_enrichment dataframe that displays the log2-transformed enrichment values
snp_enrichment <- snp_enrichment %>%
  dplyr::mutate(log2_enrichment = log2(Enrichment + 1.0))

# plot the data from the snp_enrichment dataframe using ggplot, with MAF on the x-axis and log2_enrichment on the y-axis
# geom_line - specify that the lines should be determined based on the genomic feature, not the MAF
ggplot(data = snp_enrichment, 
       mapping = aes(x = MAF, y = log2_enrichment)) + 
  geom_line(data = snp_enrichment,
            mapping = aes(color = Feature)) + 
  scale_color_colorblind(name = "Genomic Feature") +
  labs(
    title = "SNP enrichment in different genomic features",
    subtitle = "Enrichment compared across different minor allele frequencies",
    x = "Minor Allele Frequency (MAF)",
    y = "SNP Enrichment (log2)"
  )

# save the plot generated above as a pdf
ggsave(filename = "~/qbb2024-answers/week2/snp_enrichments.pdf")
