library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(ggfortify)

## Organize the data ##

# Read in the data file (in current working directory)
data = readr::read_tsv("salmon.merged.gene_counts.tsv")

# Change the row names to the name of the gene associated with that row
data = column_to_rownames(data, var = "gene_name")

# Remove the gene_id column
data = data %>% dplyr::select(-gene_id)

# Make every value in the dataframe an integer (which will allow for downstream analysis)
data = data %>% dplyr::mutate_if(is.numeric, as.integer)


## Select the data we are interested in ##

# Select only the rows/genes that have at least 100 reads
data = data[rowSums(data) > 100, ]

# Select only the "narrow region" samples
data = data %>% dplyr::select("A1_Rep1":"P2-4_Rep3")


## Create a DESeq2 model ##

# Make a metadata tibble
metadata = tibble(tissue = as.factor(c("A1", "A1", "A1", "A2-3", "A2-3", "A2-3", "Cu", "Cu", "Cu", "LFC-Fe", "LFC-Fe", "LFC-Fe", "Fe", "Fe", "Fe", "P1", "P1", "P1", "P2-4", "P2-4", "P2-4")), rep = as.factor(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)))

# Make a DESeq2 object
ddsNarrow = DESeqDataSetFromMatrix(countData = as.matrix(data), colData = metadata, design = ~tissue)

## Apply vst batch correction ##

# Look at mean by variance
meanSdPlot(assay(ddsNarrow))

# Use the vst function to batch correct
vstNarrow = vst(ddsNarrow)

# Look at mean by variance
meanSdPlot(assay(vstNarrow))


## PCA ##

# Run a PCA on the vstNarrow data
pca_data = plotPCA(vstNarrow, intgroup=c("rep","tissue"), returnData=TRUE)

# Plot the PCA
ggplot(pca_data, aes(PC1, PC2, color = tissue, shape = rep)) +
  geom_point(size=5)

# Save the plot
ggsave(filename = "PCA_plot_original.png")


## Correct the mislabeled replicates (LFC-Fe_Rep3 and Fe_Rep1) ##

# Swap the mislabeled columns
corrected_data = data %>% dplyr::relocate(Fe_Rep1, .before = `LFC-Fe_Rep3`)

# Rename the mislabeled columns
colnames(corrected_data)[12] = "LFC-Fe_Rep3"
colnames(corrected_data)[13] = "Fe_Rep1"


## Repeat the steps up to the PCA on the corrected dataframe ##

# Make a DESeq2 object
ddsNarrow_corrected = DESeqDataSetFromMatrix(countData = as.matrix(corrected_data), colData = metadata, design = ~tissue)

# Look at mean by variance
meanSdPlot(assay(ddsNarrow_corrected))

# Use the vst function to batch correct
vstNarrow_corrected = vst(ddsNarrow_corrected)

# Look at mean by variance
meanSdPlot(assay(vstNarrow_corrected))

# Run a PCA on the vstNarrow data
pca_data_corrected = plotPCA(vstNarrow_corrected, intgroup=c("rep","tissue"), returnData=TRUE)

# Plot the PCA
ggplot(pca_data_corrected, aes(PC1, PC2, color = tissue, shape = rep)) +
  geom_point(size=5)

# Save the plot
ggsave(filename = "PCA_plot_corrected.png")



## Filter genes by variance ##

# Generate a matrix from the vst-corrected data
matNarrow_corrected = as.matrix(assay(vstNarrow_corrected))

# Take the average and calculate the standard deviation for each set of replicates
combined = matNarrow_corrected[,seq(1, 21, 3)]
combined = combined + matNarrow_corrected[,seq(2, 21, 3)]
combined = combined + matNarrow_corrected[,seq(3, 21, 3)]
combined = combined / 3
sds = rowSds(combined)

# Select genes with a standard deviation greater than 1
filt = rowSds(combined) > 1
matNarrow_corrected = matNarrow_corrected[filt, ]

# Set a seed to ensure that the clustering I generate is replicable
set.seed(42)


## Apply K-means clustering to the data ##

# Apply the K-mean clustering and order the data
k = kmeans(matNarrow_corrected, centers = 12)$cluster
ordering = order(k)
k = k[ordering]
matNarrow_corrected = matNarrow_corrected[ordering, ]

# Make a saved heatmap file
png("heatmap_kmeans.png")

# Generate the heatmap and write it to the file
heatmap(matNarrow_corrected, Rowv = NA, Colv = NA, RowSideColors = RColorBrewer::brewer.pal(12, "Paired")[k])

dev.off()

# Identify and save the genes in cluster 1
cluster1 = rownames(matNarrow_corrected[k==1,])
write.table(cluster1, 'cluster1.txt', sep="\n", quote=FALSE,
            row.names=FALSE, col.names=FALSE)







