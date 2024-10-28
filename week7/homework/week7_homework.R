library(tidyverse)
library(broom)
library(DESeq2)

## Exercise 1.1 ##

# read in the data and metadata files
data <- read_delim("gtex_whole_blood_counts_downsample.txt")
metadata <- read_delim("gtex_metadata_downsample.txt")

# make the row names the gene names in the data file
data <- column_to_rownames(data, var = "GENE_NAME")

# make the row names the subject IDs in the metadata file
metadata <- column_to_rownames(metadata, var = "SUBJECT_ID")


## Exercise 1.2 ##

# make sure that the row names in the metadata file match the column names in the data file
table(colnames(data) == rownames(metadata))

# create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ SEX + DTHHRDY + AGE)


## Exercise 1.3 ##

# apply vst transformation
vsd <- vst(dds)

# apply and plot PCA and save the plots as variables
pca_sex <- plotPCA(vsd, intgroup = "SEX") + 
  labs(title = "PCA grouped by sex")
pca_dthhrdy <- plotPCA(vsd, intgroup = "DTHHRDY") + 
  labs(title = "PCA grouped by cause of death")
pca_age <- plotPCA(vsd, intgroup = "AGE") + 
  labs(title = "PCA grouped by age")

# save each of the PCA plots as PNG files
png(filename = "pca_sex.png")
plot(pca_sex)
dev.off()

png(filename = "pca_dtthrdy.png")
plot(pca_dthhrdy)
dev.off()

png(filename = "pca_age.png")
plot(pca_age)
dev.off()

# Question 1.3.3 #
# The first principle component explains 48% of the variance in the data, and the second principle component explains 7% of the variance in the data. Given those values and the attributes associated with each sample, PC1 seems to associate with cause of death, while PC2 appears to be associated with sex. This suggests that cause of death has a large effect on differential gene expression, while sex has a much smaller effect and age seems to have no measurable effect.


## Exercise 2.1 ##

# make a data frame out of the vst-transformed data
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()

# append the metadata information for each sample to the vsd data frame
vsd_df <- bind_cols(metadata, vsd_df)

# run a linear regression model on the data with the gene WASH7P as the response variable
m1 <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

# Question 2.1.1 #
# No, WASH7P does not appear to display sex-differential gene expression. The p-value calculated for males compared to females was 2.79e-1, which is far above the threshold typically needed to determine statistical significance.


# run a linear regression on the data with the gene SLC25A47 as the response variable
m2 <- lm(formula = SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()

# Question 2.1.2 #
# Yes, SLC25A47 does display sex-differential gene expression, as the p-value for males compared to females is 2.57e-2, which falls below the threshold needed to determine statistical significance. Based on the regression model, males are estimated to express SLC25A47 more highly than females, given that the slope of the SEXmale term is positive (0.518). Granted, the difference between males and females doesn't seem too extreme given the small slope.


## Exercise 2.2 ##

# perform differential analysis on the un-normalized data
dds <- DESeq(dds)


## Exercise 2.3 ##

# extract data comparing gene expression in males vs females
sex_res <- results(dds, name = "SEX_male_vs_female") %>%
  as_tibble(rownames = "GENE_NAME")

# filter for genes that display differential expression and have a FDR of 10% and arrange by padj
sex_res_filt <- sex_res %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1) %>%
  arrange(padj)

# return the number of genes that fit the filter requirements above
sex_res_filt %>% nrow()

# Question 2.3.2 #
# There are 262 genes that exhibit significant differential expression between males and females.


# read in the gene locations file
fb <- read_delim("gene_locations.txt")

# append the chromosome information for each gene to the sex_res tibble
sex_res <- left_join(sex_res, fb, by = "GENE_NAME") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# Question 2.3.3 #
# The top 11 genes (and 17 out of the top 20 genes) that display differential expression between males and females are all located on the Y chromosome. Not surprisingly, since only males have a Y chromosome, these genes are expressed to a greater extent in males than in females. 3 out of the top 20 genes are located on the X chromosome. Since females have two X chromosomes, it makes sense that these genes display lower expression in males than in females. Overall, it is evident that there are more male-upregulated genes near the top of the list, which is probably an artifact of the fact that males still have an X chromosome, which makes differential expression for genes on the X chromosome more difficult to achieve.


# return the information for the gene WASH7P and the gene SLC25A47
sex_res %>% filter(GENE_NAME == "WASH7P" | GENE_NAME == "SLC25A47")

# Question 2.3.4 #
# Yes, the results are largely consistent between the simple linear regression analysis run earlier and the DESeq differential expression analysis. Specifically, in both tests, WASH7P was not shown to display statistically significant differential expression between males and females, whereas SLC25A47 was.


## Exercise 2.4 ##

# extract data comparing gene expression in patients with different causes of death
death_res <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes") %>%
  as_tibble(rownames = "GENE_NAME")

# filter for genes that display differential expression and have a FDR of 10% and arrange by padj
death_res <- death_res %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.1) %>%
  arrange(padj)

# return the number of genes that fit the filter requirements above
death_res %>% nrow()

# Question 2.4.1 #
# 16,069 genes are differentially expressed based on death classification.

# Question 2.4.2 #
# The fact that more genes are differentially expressed based on cause of death rather than sex makes sense given the qualitative analysis I performed in exericse 1.3 (question 1.3.3). According to the PC analysis and based on my interpretation of the plots I generated, cause of death accounted for about 48% of the variance in gene expression across subjects, while sex only explained about 7% of the variance in gene expression across the population.


## Exercise 3 ##

# make a volcano plot to highlight genes that are differentially expressed based on sex
ggplot(data = sex_res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = (abs(log2FoldChange) > 1 & padj < 1e-1))) +
  xlim(-7, 10) +
  ylim(0, 300) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray") +
  geom_vline(xintercept = 1, linetype = 2, color = "gray") +
  geom_vline(xintercept = -1, linetype = 2, color = "gray") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkgray", "darkred")) +
  labs(title = "Differential gene expression in whole blood samples based on sex",
       y = expression(-log[10]("q-value")), 
       x = expression(log[2]("fold change")))

# save the volcano plot
ggsave(filename = "deseq_analysis_volcano.png")



