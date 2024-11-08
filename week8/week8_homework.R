# Exercise 1 - Load Data #

# load in packages and libraries
# BiocManager::install("zellkonverter")
library("zellkonverter")
library("scater")
library("scran")
library("scuttle")
library("ggplot2")
library("scales")

# read in the data file and normalize
gut <- readH5AD("~/qbb2024-answers/week8/data/v2_fca_biohub_gut_10x_raw.h5ad")
assayNames(gut) <- "counts"
gut <- logNormCounts(gut)


## Question 1 ##

# view gut
gut

# There are 13,407 genes quantitated
# There are 11,788 cells in the dataset
# There are PCA, tSNE, and UMAP dimension reduction datasets



## Question 2 ##

# assess colnames in gut
as.data.frame(colData(gut))
colnames(as.data.frame(colData(gut)))

# plot the UMAP for gut, colored based on broad annotation
plotReducedDim(gut, "X_umap", color_by = "broad_annotation")

# There are 39 columns
# The most interesting columns to me are n_counts, n_genes, and sex, because it seems like those data points would be the most important to consider for differential expression analysis.




# Exercise 2 - Explore Data #

## Question 3 ##

# sum expression of each gene, plot the distribution of the sums, and determine which genes are most highly expressed
genecounts <- rowSums(assay(gut))
summary(genecounts)
hist(genecounts)
head(sort(genecounts, decreasing = TRUE))

# Mean gene count is 3185 and median gene count is 254. The lower median suggests that most of the samples have low gene counts, but the higher mean suggests that some of the samples have incredibly high gene counts.
# The three genes with highest expression are Hsromega (lncRNA), CR45845 (pre-rRNA), and roX1 (lncRNA). These are all non protein coding RNAs



## Question 4a ##

# sum the total expression for all genes in each cell and plot the distribution of sums
cellcounts <- colSums(assay(gut))
summary(cellcounts)
hist(cellcounts)

# Mean number of counts per cell is 3622.
# Cells with high total counts must be specialized cell types in the gut that express a high volume of protein. This could represent digestive cells that release high amounts of enzymes for digestion.



## Question 4b ##

# sum the number of genes detected above the threshold of 0 for each cell and plot the distribution of sums
celldetected <- colSums(assay(gut) > 0)
summary(celldetected)
hist(celldetected)

# Mean number of genes detected per cell is 1059.
# This mean represents 0.078, or 7.8% (1059 / 13,407) of the total number of genes in the dataset.



## Question 5 ##

# sum expression of mitochondrial genes and add the data to the original single cell experiment object gut
rownames(gut)
mito <- grep(rownames(gut), pattern = "^mt:", value = TRUE)
df <- perCellQCMetrics(x = gut, subsets = list(Mito=mito))
df <- as.data.frame(df)
summary(df)
colData(gut) <- cbind( colData(gut), df )

# plot the percent of reads that map to mitochondrial genes based on the broad annotations
plotColData(gut, y = "subsets_Mito_percent", x = "broad_annotation") + 
  theme( axis.text.x=element_text( angle=90 ) )

# save the plot as a png
png(filename="~/qbb2024-answers/week8/broad_annotations_vs_percent_mito.png")
plotColData(gut, y = "subsets_Mito_percent", x = "broad_annotation") + 
  theme( axis.text.x=element_text( angle=90 ) )
dev.off()

# Epithelial cells seem to have the most cells with a higher percent of mitochondrial reads. This may make sense as epithelial cells likely need to grow and divide relatively rapidly in the gut, which requires ATP and, therefore, requires more mitochondrial proteins to be expressed. A similar reasoning could explain why somatic precursor cells seem to have the highest average percent mitochondrial reads.




# Exercise 3 - Identify Marker Genes #

## Question 6a ##

# subset the gut object to include only cells that are annotated as epithelial cells. Plot the UMAP for these epithelial cells and color based on specific cell type
coi <- colData(gut)$broad_annotation == "epithelial cell"
epi <- gut[,coi]
plotReducedDim(epi, "X_umap", color_by = "annotation")

# save the UMAP as a png
png(filename="~/qbb2024-answers/week8/epithelial_UMAP.png")
plotReducedDim(epi, "X_umap", color_by = "annotation")
dev.off()



## Question 6b ##

# identify the top marker genes in the anterior midgut
marker.info <- scoreMarkers( epi, colData(epi)$annotation )
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])

# The top marker genes in the anterior midgut of the adult fly are Mal-A6, Men-b, vnd, betaTry, Mal-A1, and Nhe2. Two of the genes (Mal-A6 and Mal-A1) encode glucosidases, which are involved in carbohydrate metabolism. One of the genes (betaTry) encodes trypsin, which is involved in protein metabolism.

# plot expression of the top marker gene in each specific epithelial cell type
plotExpression(epi, "Mal-A6", x = "annotation") +
  theme( axis.text.x=element_text( angl=90, vjust = 0.5) ) + 
  scale_x_discrete(labels = label_wrap(25))

# save the plot as a png
png(filename = "~/qbb2024-answers/week8/epithelial_expression_plot.png")
plotExpression(epi, "Mal-A6", x = "annotation") +
  theme( axis.text.x=element_text( angl=90, vjust = 0.5) ) + 
  scale_x_discrete(labels = label_wrap(25))
dev.off()



## Question 7 ##

# subset the gut object to include only cells that are annotated as somatic precursor cells and identify the top marker genes in the intestinal stem cells
coi <- colData(gut)$broad_annotation == "somatic precursor cell"
spc <- gut[,coi]
marker.info <- scoreMarkers( spc, colData(spc)$annotation )
chosen <- marker.info[["intestinal stem cell"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing=TRUE),]
head(ordered[,1:4])

# Marker genes for intestinal cells include hdc, kek5, N, zfh2, Tet, and Dl

# plot expression of the top six marker genes for each somatic precursor cell type
goi <- rownames(ordered)[1:6]
plotExpression(spc, goi, x = "annotation") +
  theme( axis.text.x=element_text( angl=45, hjust = 1) ) + 
  scale_x_discrete(labels = label_wrap(25))

# save the plot as a png
png(filename = "~/qbb2024-answers/week8/marker_genes_precursor_cells.png")
plotExpression(spc, goi, x = "annotation") +
  theme( axis.text.x=element_text( angl=45, hjust = 1) ) + 
  scale_x_discrete(labels = label_wrap(25))
dev.off()

# Enteroblasts and intestinal stem cells seem to have more similar expression of these top 6 marker genes.
# Given the large number of intestinal stem cells that highly express it, Dl seems to be the most specific marker for intestinal stem cells.


