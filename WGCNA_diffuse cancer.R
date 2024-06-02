#CSB project
# Load required libraries
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(flashClust)
library(dplyr)  

allowWGCNAThreads() 

# 1. Fetch Data ---------------------


# Read the counts data into R
data <- read.csv("~/Desktop/GSE113255_counts 2.csv", header = TRUE, row.names = 1)

# Check the structure of the data
str(data)

# get metadata
geo_id <- "GSE113255"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(2,8:13,43,44)]

# Define the sample IDs from your gene expression data
sample_ids <- c("GSM3101081", "GSM3101062", "GSM3101063", "GSM3101064", "GSM3101085", "GSM3101086", "GSM3101067", "GSM3101068", "GSM3101087", "GSM3101070", "GSM3101071", "GSM3101088", "GSM3101072", "GSM3101074", "GSM3101075", "GSM3101076", "GSM3101077", "GSM3101078", "GSM3101092", "GSM3101080", "GSM3101103", "GSM3101104", "GSM3101105", "GSM3101090", "GSM3101093", "GSM3101106", "GSM3101118", "GSM3101120", "GSM3101129", "GSM3101134")

# Subset the metadata to keep only the rows with sample IDs
phenoData <- phenoData[grepl(paste(sample_ids, collapse = "|"), rownames(phenoData)), ]

# Change column name
colnames(phenoData)[colnames(phenoData) == "cancer type:ch1"] <- "cancer type"

# Check the updated metadata
head(phenoData)


# 2. QC - outlier detection 
# detect outlier genes
# goodSamplesGenes is from wcgna package 

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

# dont include this pca as its messy (rashu)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
data.subset <- data
colData <- phenoData 
names(colData)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model


## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) 


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
assay(dds_norm) %>% 
  head()

# we transposed the normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
power

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 12
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 10000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Calculate adjacency and TOM
adjacency <- adjacency(norm.counts, power = 12)
TOMadj <- TOMsimilarity(adjacency)
dissTOMadj <- 1 - TOMadj


# Clustering using TOM
# Call the hierarchical clustering function 
hclustGeneTree <- hclust(as.dist(dissTOMadj), method = "average")

# Plot the resulting clustering tree (dendogram)
sizeGrWindow(12, 9)
plot(hclustGeneTree, xlab = "", sub = "", 
     main = "Gene Clustering on TOM-based disssimilarity", 
     labels = FALSE, hang = 0.04)

# Make the modules larger, so set the minimum higher
minModuleSize <- 30

# Module ID using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = hclustGeneTree, 
                             distM = dissTOMadj,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

table(dynamicMods)

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(hclustGeneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
dynamic_MEList <- moduleEigengenes(norm.counts, colors = dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes

# Calculate dissimilarity of module eigengenes
dynamic_MEDiss <- 1-cor(dynamic_MEs)
dynamic_METree <- hclust(as.dist(dynamic_MEDiss))
# Plot the hclust
sizeGrWindow(7,6)
plot(dynamic_METree, main = "Dynamic Clustering of module eigengenes",
     xlab = "", sub = "")

######################## MERGE SIMILAR MODULES
dynamic_MEDissThres <- 0.50

# Plot the cut line
#abline(h = dynamic_MEDissThres, col = "red")

# Call an automatic merging function
merge_dynamic_MEDs <- mergeCloseModules(norm.counts, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)
# The Merged Colors
dynamic_mergedColors <- merge_dynamic_MEDs$colors

# Eigen genes of the new merged modules
mergedMEs <- merge_dynamic_MEDs$newMEs
mergedMEs

table(dynamic_mergedColors)

sizeGrWindow(12,9)
plotDendroAndColors(hclustGeneTree, cbind(dynamicColors, dynamic_mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


library(dplyr)

# Create a new data frame called traits from phenoData
traits <- phenoData

# Binarize characteristics_ch1.1 column
traits$char_ch1.1_binary <- ifelse(traits$characteristics_ch1.1 == "cancer type: diffuse type", 1, 
                                   ifelse(traits$characteristics_ch1.1 == "cancer type: n/a", 0, NA))

# Binarize characteristics_ch1 column
traits$char_ch1_binary <- ifelse(traits$characteristics_ch1 == "tissue: gastric cancer", 0,
                                 ifelse(traits$characteristics_ch1 == "tissue: normal intestinal mucosae", 1, NA))

# Print the updated traits
print(traits)
# Remove the column "gender:ch1" from the traits dataframe
traits <- select(traits, -`gender:ch1`)

# Print the updated traits dataframe
print(traits)
# Rename columns
traits <- traits %>%
  rename(Diffuse_GC = char_ch1.1_binary,
         Normal = char_ch1_binary)
# Select columns 9 and 10
traits <- select(traits, 9:10)
# Print the updated traits dataframe
print(traits)

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Calculate correlations
module.trait.corr <- cor(mergedMEs, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(mergedMEs, traits, by = 'row.names')

head(heatmap.data)
names(heatmap.data)


CorLevelPlot(heatmap.data,
            
             x = names(heatmap.data)[26:27],
             y = names(heatmap.data)[2:25],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(mergedMEs, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]

# Calculate the gene significance and associated p-values for Diffuse_GC
gene.signf.corr_diffuse <- cor(norm.counts, traits$Diffuse_GC, use = 'p')
gene.signf.corr.pvals_diffuse <- corPvalueStudent(gene.signf.corr_diffuse, nSamples)

# Calculate the gene significance and associated p-values for Normal
gene.signf.corr_normal <- cor(norm.counts, traits$Normal, use = 'p')
gene.signf.corr.pvals_normal <- corPvalueStudent(gene.signf.corr_normal, nSamples)

# Combine the p-values into a dataframe
pvals_combined <- data.frame(Diffuse_GC = gene.signf.corr.pvals_diffuse, 
                             Normal = gene.signf.corr.pvals_normal)

# Print the top rows of the combined p-values dataframe
print(head(pvals_combined))

# Print the top 25 rows of the combined p-values dataframe
print(head(pvals_combined, 25))

# Extract gene IDs from the row names of the top 25 rows
top_25_gene_ids <- rownames(head(pvals_combined, 25))

# Print the extracted gene IDs
print(top_25_gene_ids)

library(org.Hs.eg.db)
# Your top 25 gene IDs
top_25_gene_ids <- c("653635", "729737", "102723897", "100132287", "113219467", "100133331", "100288069", "105378580",
                     "643837", "148398", "26155", "339451", "84069", "57801", "9636", "375790", "54991", "8784",
                     "7293", "51150", "126792", "118424", "116983", "126789", "54973")

# Convert gene IDs to gene symbols
top_25_gene_symbols <- mapIds(org.Hs.eg.db, keys = top_25_gene_ids, column = "SYMBOL", keytype = "ENTREZID")

print(top_25_gene_symbols)

