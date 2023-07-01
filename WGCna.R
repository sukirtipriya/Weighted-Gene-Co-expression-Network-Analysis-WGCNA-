# Script to perform WGCNA

library(WGCNA)
library(DESeq2)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()

# 1.Fetch Data

data <- read.delim("GSE152418_p20047_Study1_RawCounts.txt", header=T)

# Get metadata

geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1,2,46:50)]

# Prepare Data

data[1:10,1:10]

data <- data %>%
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>%
  mutate(samples = gsub('\\_','-',samples)) %>%
  inner_join(., phenoData, by = c('samples' = 'title')) %>%
  select(1,3,4) %>%
  spread(key = "geo_accession", value = "counts") %>%
  column_to_rownames(var = "ENSEMBLID")

head(data)
  

# 2.QC - outliers detection
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allok

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected as outlier

data <- data[gsg$goodGenes ==TRUE,]

# detect outlier samples - hierarchial clustering - method 1

htree <- hclust(dist(t(data)), method = "average")
plot(htree)

#pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca,dat)

ggplot(pca.dat, aes(PC1, PC2))+
  geom_point()+
  geom_text(label = rownames(pca.dat))+
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2], '%'))
       
# exclude outlier samples

samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# 3.Normalization

# create a deseq2 dataset
# exclude outlier samples

colData <- phenoData %>%
  filter(!row.names(.) %in% samples.to.be.excluded)

# fixing column names in colData

names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column name identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

#create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design =~ 1)
 
dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]

# perform variance stablization
dds_none <- vst(dds75)

#get normalized counts
norm.count <- assay(dds_norm) %>%
   t()

# 4. Network construction
 
 # choose a set of soft-thersholding powers
 power <- c(c(1:10), seq(from = 12, to = 50, by =2))
 
 #call the network topology analysis function
 sft <- pickSoftThreshold(norm.counts,
                    powerVector = power,
                    networkType = "signed",
                    verbose = 5)


sft.data <- sft$fitIndices

# visualization t pick power

a1 <- ggplot(sft.data, aes(power, SFT.R.sq, label = power))+
          geom_point() +
          geom_text(nudge_y = 0.1)+
          geom_hline(yintercept = 0.8, color = 'red') +
          labs(x = 'power', y = 'scale free topology model fit, signed R^2')+
          theme_classic()
          

 a2 <- ggplot(sft.data, aes(power, mean.k., label = power))+
          geom_point() +
          geom_text(nudge_y = 0.1)+
          labs(x = 'power', y ='Mean connectivity')+
          theme_classic()

grid.arrange(a1,a2,norm = 2)


# convert matrix to numeric

norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize

bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = 'signed',
                          power = 'soft_power',
                          mergeCutHeight = 0.25,
                          numericlabels = FALSE,
                          ramdomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

# 5. Module Eigengenes

module_eigengenes <- bwnet$MEs

head(module_eigengenes)

# get number of genes for each module

table(bwnet$colors)

# plot the dendogram and the module colors before and after merging underneath

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged","merged"),
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHong = 0.05)
                    
# 6. Realate module to traits

#create traits file - binorize categorical variables

traits <- colData %>% 
               mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1,0)) %>%
               select(7)

# binarize categorical variables

colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", 
                          "Moderate", "Severe"))

severity.out <- binarizeCategoricalColumns(colData$severity,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)
                           
                           
traits <- cbind(traits, severity.out)


# Define number of genes and samples

nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualization module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>%
                       column_to_rownames(var = 'Row.names')


corLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22],
             y = names(heatmap.data)[1:17],
             col = c("blue1","skyblue","white","pink","red"))
             
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
          filter('bwnet$colors' == 'turquoise') %>%
          rownames()

# Interamodular analysis: Identify driver genes

module.membership.measure <- cor(module_eigengenes, norm.counts, use ='p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]

# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use ='p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>%
                as.data.frame() %>%
                arrange(V1) %>%
                head(25)
                

