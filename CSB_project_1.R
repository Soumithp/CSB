
install.packages("remotes")
remotes::install_github("kevinblighe/CorLevelPlot")



library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

data<- read.csv("C:/Users/sowmi/Desktop/4th Semester/Systems_biology/project_apr26/GSE202203_RawCounts_gene_3207.tsv", header = T, sep = "\t")

data[,1]

#gettingmetadata
geo_id <- "GSE202203"

gse<- getGEO(geo_id, GSEMatrix = T)


phenotype_data1 <- pData(phenoData(gse[[1]]))

phenotype_data2 <- pData(phenoData(gse[[2]]))



combined_phenotypes <- rbind(phenotype_data1, phenotype_data2 )

combined_phenotypes<- combined_phenotypes[, c(1,2,65,66,67,69,72,82,87)]



row.names(combined_phenotypes) <- combined_phenotypes[, "title"]

combined_phenotypes<-combined_phenotypes[,-1]




#preparing the data 
colnames(data)

row.names(data) <- data[, "X"]


data<-data[,-1]


colnames(data)


#QC-detecting outliers 
gsg<- goodSamplesGenes(t(data))

data <- data[gsg$goodGenes == TRUE,]

h_tree<-hclust(dist(data), method = "average")
plot(h_tree)


#pcamethod for outliers
pca<- prcomp(t(data))

pca_data<- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)




ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


samples.to.be.excluded <- c("S000176", "S000625", "S000985", "S000795", "S001071", "S001163", "S001395", "S001810")

data_excluded_OL <- data[, !(colnames(data) %in% samples.to.be.excluded)]


#normalization of the data

# exclude outlier samples


colData <- combined_phenotypes %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)

# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))



# making the rownames and column names identical
all(rownames(colData) %in% colnames(data_excluded_OL))


length(setdiff(rownames(colData), colnames(data_excluded_OL))) == 0
all(sort(rownames(colData)) == sort(colnames(data_excluded_OL)))


all(rownames(colData) == colnames(data_excluded_OL))

idx <- match(rownames(colData), colnames(data_excluded_OL))
countData <- data_excluded_OL[, idx]


# Match row names and reorder countData

#data_excluded_OL<- t(data_excluded_OL)

countData<- round(countData)


# create a deseq2 dataset
# create dds
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ 1) 

## remove all genes with counts < 15 in more than 75% of samples
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 3200,]
nrow(dds75) # 8780 genes


# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
norm.counts <- assay(dds_norm) %>%  t()

class(norm.counts)
str(norm.counts)

# 4. Network Construction
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#norm.counts <- norm.counts[complete.cases(norm.counts),]

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

#visualization for selecting threshold

Rsquareplot <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

#for mean
meanplot <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(Rsquareplot, meanplot, nrow = 2)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 5
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
cor <- temp_cor

# 5. Module Eigengenes-
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)

# get number of genes for each module
table(bwnet$colors)


# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# 6A. Relate modules to traits-
# module trait associations


# create traits file - binarize categorical variables
traits <- colData %>% 
  select(chemo_treated)

library(dplyr)

# Store the row names of the dataframe
row_names_original <- rownames(traits)


# Convert the "chemo_treated" column to numeric
traits$chemo_treated <- as.numeric(traits$chemo_treated)

# Replace the NA values in the "chemo_treated" column with 0
traits$chemo_treated <- ifelse(is.na(traits$chemo_treated), 0, traits$chemo_treated)


# Assign the previously stored row names to the dataframe
rownames(traits) <- row_names_original

#traits$chemo_untreated=0


colData$age <- ifelse(colData$age_at_diagnosis >= 24 & colData$age_at_diagnosis <= 50, 1,
                      ifelse(colData$age_at_diagnosis > 50 & colData$age_at_diagnosis <= 96, 2, NA))

colData <- colData[!is.na(colData$age), ]
age_data <- data.frame(sample_names = rownames(colData), age = colData$age)


# Save the values of the "sample" column to a vector
sample_names <- age_data$sample_names

# Remove the "sample" column from the "age_data" dataframe
age_data <- age_data[-which(names(age_data) == "sample_names")]

# Assign the saved "sample" column values as row names to the "age_data" dataframe
rownames(age_data) <- sample_names

traits <- cbind(traits, age_data)


treated = as.data.frame(traits$chemo_treated);



MET = orderMEs(cbind(module_eigengenes, treated))

# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
module = "red"
modNames = substring(names(module_eigengenes), 3)

column = match(module, modNames);


#
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(module.membership.measure[module.gene.mapping, column]),
                   
                   abs(gene.signf.corr[module.gene.mapping, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)



module.trait.corr <- cor(module_eigengenes, traits, use = "p")

module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[36:37],
             y = names(heatmap.data)[1:35],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# Display the correlation values within a heatmap
labeledHeatmap(Matrix = module.trait.corr, 
               xLabels = colnames(traits),
               yLabels = names(module_eigengenes), 
               ySymbols = names(module_eigengenes), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait Relationships"))



module.gene.mapping <- as.data.frame(bwnet$colors)

module.gene.mapping %>% 
  filter(`bwnet$colors` == 'red') %>% 
  rownames()


#6B. Intramodular analysis: Identifying driver genes ---------------
  
  
  
  # Calculate the module membership and the associated p-values
  
  # The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
  # This quantifies the similarity of all genes on the array to every module.
  
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]



# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$chemo_treated, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

selected_rows <- heatmap.data[ , "MEred"]
