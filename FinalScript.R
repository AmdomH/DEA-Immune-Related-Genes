
#Load library
library(tidyverse) 
library(DESeq2)
library(readxl)
library(pheatmap)
library(psych) # for descriptive statistics
library(EnhancedVolcano) # for enhanced volcano plot
library(ggthemes)
if(!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("GO.db", quietly = TRUE))  BiocManager::install("GO.db", update = TRUE, ask = FALSE)
library(WGCNA)
library(GEOquery)
library(ComplexHeatmap) # required for ComplexHeatmap
library(EnhancedVolcano)# for nicer volcano plots
library(UpSetR) #bar plots
library(gridExtra)
library(CorLevelPlot)

#Set working directory
setwd("C:/Data")
#Load count data and Meta data
# Load count_data from .Rdata file
load("count_data.Rdata")
load("meta.Rdata")
load("gene.Rdata")
load("immune_genes.Rdata")

###################################
# Chapter 01: Descriptive Analysis
####################################
ggplot(meta, aes(x=diagnosis, fill=gender)) +
  geom_bar(position="dodge") +
  labs(x="Diagnosis", y="Count", fill="Gender") +
  theme_solarized() 

meta%>%
  group_by(gender)%>%
  summarise(perc= n()/nrow(meta) * 100)

#create table
table(meta$diagnosis, meta$age_cat)

#Box plot
ggplot(meta, aes(x=diagnosis, y=age, fill=diagnosis)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, color = "black")+
  ylab("Age in Years")+
  labs( x="Diagnosis", y="Age", fill="Diagnosis") +
  scale_y_continuous(breaks= c(0:20))+
  theme_economist_white()

################################################
## Chapter 02 Exploratory Analysis
###############################################
dds<- DESeqDataSetFromMatrix(countData = count_data,
                             colData = meta,
                             design = ~ diagnosis)
dds_normal<- estimateSizeFactors(dds)


sizeFactors(dds_normal)[1:5] #check top 5 samples
count_normal<- counts(dds_normal, normalized = TRUE) #Extrating Normalized counts
count_normal[1:5,1:5]


vsd <- vst(dds_normal, blind = TRUE)#scaling the data
vsd_mat<- assay(vsd)#Extract the vst matrix
vsd_mat[1:3, 1:3]

#subset the immune related genes
by <- join_by(gene_name==GeneName)
immune_genes<- left_join(immune_genes, gene,by)
immune_genes%>%head()

vsd_immune_mat <- vsd_mat[rownames(vsd_mat) %in% immune_genes$ID.x, ]
dim(vsd_immune_mat)
vsd_immune_cor<- cor(vsd_immune_mat) #compute correlation
vsd_immune_cor[1:5,1:5]

set.seed(123)
# Create the Heatmap object
p<-Heatmap(vsd_immune_cor, 
             name = "cor", 
             show_column_names = FALSE,  # hide column labels
             show_row_names = FALSE,  # hide row labels
             cluster_rows = TRUE,  # cluster rows
             cluster_columns = TRUE,  # cluster columns
             top_annotation = HeatmapAnnotation(diagnosis = meta$diagnosis),
             column_title_rot = 90
)

draw(p, heatmap_legend_side = "right")# Draw the heatmap
#Plot PCA
pcaData <- plotPCA(vsd, intgroup = "diagnosis", returnData = TRUE)
cbbPalette <- c("#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
pcaData%>%
  ggplot(aes(PC1, PC2, color=diagnosis))+
  geom_point(size = 2, alpha=1)+
  scale_color_manual(values =cbbPalette)+
  theme_economist_white()

################################################
# Chapter 03 Differential analysis
#################################################
dds <- DESeq(dds)# Perform differential expression analysis     
norm_counts <- counts(dds, normalized = TRUE)
norm_counts[1:5, 1:5]

plotDispEsts(dds)
#Checking the distributin count data
ggplot(count_data) +
  geom_histogram(aes(x = log2(sample_63+1)), stat = "bin", bins = 30) +
  xlab("expression counts") +
  ylab("Number of genes") +
  theme_economist_white()

#Plot mean-variance
mean_counts <- apply(count_data[,1:147], 1, mean)
variance_counts <- apply(count_data[,1:147], 1, var)
df <- data.frame(mean_counts, variance_counts)
p<-ggplot(df) +
  geom_point(aes(x=log2(mean_counts), y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")+
  theme_economist_white()
p

#create contrasts
diagnosis <- c('Craniopharyngioma', 'ATRT', 'Ependymoma', 'Glioblastoma', 'Glioma', 'Medulloblastoma')
contrasts <- lapply(setdiff(diagnosis, "Craniopharyngioma"), function(x) c("diagnosis", x, "Craniopharyngioma"))
contrasts
results.list <- lapply(contrasts, function(con) results(dds, contrast = con, alpha = 0.05))

#### ATRT vs Craniopharyngioma
ATRT_res<-results.list[[1]]
ATRT_res
summary(ATRT_res)
#### Ependymoma vs Craniopharyngioma
Ependymoma_res<- results.list[[2]]
Ependymoma_res
summary(Ependymoma_res)
#### Glioblastoma vs craniopharungioma
Glioblastoma_res<- results.list[[3]]
Glioblastoma_res
summary(Glioblastoma_res)
#### Glioma vs Craniopharyngioma
Glioma_res<- results.list[[4]]
Glioma_res
summary(Glioma_res)
#### Medulloblastoma vs Craniopharyngioma
Medulloblastoma_res<- results.list[[5]]
Medulloblastoma_res
summary(Medulloblastoma_res)

## Add gene names
#idx <- match( rownames(res), gene$ID )
ATRT_res$geneName<- gene$GeneName[match(rownames(ATRT_res), gene$ID )]
Ependymoma_res$geneName<- gene$GeneName[match(rownames(Ependymoma_res), gene$ID )]
Glioblastoma_res$geneName<- gene$GeneName[match(rownames(Glioblastoma_res), gene$ID )]
Glioma_res $geneName<- gene$GeneName[match(rownames(Glioma_res ), gene$ID )]
Medulloblastoma_res$geneName<- gene$GeneName[match(rownames(Medulloblastoma_res), gene$ID )]
head(ATRT_res)
#Filter immune related genes
idx <- which(ATRT_res$geneName %in% immune_genes$gene_name)#create index
ATRT_immune_res <- ATRT_res[idx,]
Ependymoma_immune_res<- Ependymoma_res[idx,]
Glioblastoma_immune_res<- Glioblastoma_res[idx,]
Glioma_immune_res<- Glioma_res[idx,]       
Medulloblastoma_immune_res<- Medulloblastoma_res[idx,]
head(Medulloblastoma_immune_res)


## Summary results
summary(ATRT_immune_res );
head(ATRT_immune_res)
summary(Ependymoma_immune_res);
head(Ependymoma_immune_res)
summary(Glioblastoma_immune_res)
head(Glioblastoma_immune_res)
summary(Glioma_immune_res);
head(Glioma_immune_res)
summary(Medulloblastoma_immune_res);
Medulloblastoma_immune_res

################################################
# Chapter 04 Filtering DE genes
################################################
# ATRT
top_upregulated_ATRT <- ATRT_immune_res %>%as.data.frame()%>%
  arrange(desc(log2FoldChange)) %>%
  head(5)
top_downregulated_ATRT <- ATRT_immune_res %>%as.data.frame()%>%
  arrange(log2FoldChange) %>%
  head(5)
# Ependymoma
top_upregulated_Ependymoma <- Ependymoma_immune_res %>%as.data.frame()%>%
  arrange(desc(log2FoldChange)) %>%
  head(5)
top_downregulated_Ependymoma <- Ependymoma_immune_res %>%as.data.frame()%>%
  arrange(log2FoldChange) %>%
  head(15)
# Glioblastoma
top_upregulated_Glioblastoma <- Glioblastoma_immune_res %>%as.data.frame()%>%
  arrange(desc(log2FoldChange)) %>%
  head(5)
top_downregulated_Glioblastoma <- Glioblastoma_immune_res %>%as.data.frame()%>%
  arrange(log2FoldChange) %>%
  head(5)
# Glioma
top_upregulated_Glioma <- Glioma_immune_res %>%as.data.frame()%>%
  arrange(desc(log2FoldChange)) %>%
  head(5)
top_upregulated_Glioma
top_downregulated_Glioma <- Glioma_immune_res %>%as.data.frame()%>%
  arrange(log2FoldChange) %>%
  head(5)
top_downregulated_Glioma
# Medulloblastoma
top_upregulated_Medulloblastoma <- Medulloblastoma_immune_res %>%as.data.frame()%>%
  arrange(desc(log2FoldChange)) %>%
  head(5)
top_downregulated_Medulloblastoma <- Medulloblastoma_immune_res %>%as.data.frame()%>%
  arrange(log2FoldChange) %>%
  head(5)
#################################################
# Chapter 05 Visualizations
################################################
#Volcano plot
p1<- EnhancedVolcano(ATRT_immune_res,
                     lab = ATRT_immune_res$geneName,
                     x = 'log2FoldChange',
                     y = 'padj',
                     pCutoff = 0.05,
                     FCcutoff = 2,
                     cutoffLineType = 'dotted',
                     labSize = 3,
                     title=NULL,
                     xlim = c(-30, 30),
                     selectLab = c(top_upregulated_ATRT$geneName,top_downregulated_ATRT$geneName))
p1 + theme_economist_white()
#MA Plot
plotMA(ATRT_immune_res, ylim=c(-2,2))
#Adjusted P-values plots
ggplot(as.data.frame(ATRT_immune_res), aes(padj)) +
  geom_histogram(color = "black", fill = "#0072B2", bins = 30) +
  geom_vline(aes(xintercept = 0.05), color = "red",size = 1) +
  geom_hline(aes(yintercept = 20), color = "red",size = 1) +
  geom_segment(aes(x = 0.2, y = 100, xend = 0.001, yend = 10), arrow = arrow(length = unit(0.1, "cm")), size=1.5, color="orange") +
  geom_text(aes(x = 0.2, y = 100, label = "FDR"), hjust = -0.1, vjust = 1.5, size = 5, color = "black")
labs(x = "padj", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_minimal() 

# Get DEGs for each result
# Order results based on padj
ATRT_immune_res_filter <- subset(ATRT_immune_res, padj < 0.05 & abs(log2FoldChange) > 2)
Ependymoma_immune_res_filter <- subset(Ependymoma_immune_res, padj < 0.05 & abs(log2FoldChange) > 2)
Glioblastoma_immune_res_filter <- subset(Glioblastoma_immune_res, padj < 0.05 & abs(log2FoldChange) > 2)
Glioma_immune_res_filter <- subset(Glioma_immune_res, padj < 0.05 & abs(log2FoldChange) > 2)
Medulloblastoma_immune_res_filter <- subset(Medulloblastoma_immune_res, padj < 0.05 & abs(log2FoldChange) > 2)
ATRT_immune_res_filter <- ATRT_immune_res_filter[order(ATRT_immune_res_filter$padj),]
Ependymoma_immune_res_filter <- Ependymoma_immune_res_filter[order(Ependymoma_immune_res_filter$padj),]
Glioblastoma_immune_res_filter <- Glioblastoma_immune_res_filter[order(Glioblastoma_immune_res_filter$padj),]
Glioma_immune_res_filter <- Glioma_immune_res_filter[order(Glioma_immune_res_filter$padj),]
Medulloblastoma_immune_res_filter <- Medulloblastoma_immune_res_filter[order(Medulloblastoma_immune_res_filter$padj),]
ATRT_immune_res_filter 
ATRT_DEGs <- ATRT_immune_res_filter$geneName
Ependymoma_DEGs <- Ependymoma_immune_res_filter$geneName
Glioblastoma_DEGs <-Glioblastoma_immune_res_filter$geneName
Glioma_DEGs <- Glioma_immune_res_filter$geneName
Medulloblastoma_DEGs <- Medulloblastoma_immune_res_filter$geneName

# Create a list of DEGs
degs <- list(ATRT = ATRT_DEGs, Ependymoma = Ependymoma_DEGs, Glioblastoma = Glioblastoma_DEGs,
             Glioma = Glioma_DEGs, Medulloblastoma = Medulloblastoma_DEGs)
# Make the UpSet plot
upset(fromList(degs), order.by = "freq", set_size.show=FALSE, matrix.color = "blue", text.scale = 1,
      sets.bar.color =c("maroon", "blue", "green", "purple", "orange")
)

#################################################
# Chapter 06 WGCNA
################################################
count_data_filter <- count_data %>% 
  rownames_to_column("ENSEMBL") %>%  
  filter(ENSEMBL %in% rownames(ATRT_immune_res)) %>%  
  column_to_rownames("ENSEMBL")  
count_data_filter[1:3,1:3]

gsg <- goodSamplesGenes(t(count_data_filter)) #Explore good samples
summary(gsg)
gsg$allOK

table(gsg$goodGenes);table(gsg$goodSamples)

data <- count_data_filter[gsg$goodGenes == TRUE,]# remove genes that are detectd as outliers
#PCA
pca <- prcomp(t(data))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
#create this as data.frame
pca.dat <- as.data.frame(pca.dat)
ggplot(pca.dat, aes(PC1, PC2, color= )) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

samples.to.be.excluded <- c("sample_81")#Outliers
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
colData <- meta %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)

# create dds
dds2 <- DESeqDataSetFromMatrix(countData = data.subset,
                               colData = colData,
                               design = ~ 1) # not specifying model
## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ
dds75 <- dds2[rowSums(counts(dds2) >= 15) >= 24,]
nrow(dds75) #

# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


## 4 Network Construction: The main network
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
#### Visualize 
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


##4A convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)
soft_power <- 4
temp_cor <- cor
cor <- WGCNA::cor
## Co-expression
# Calculate the co-expression (correlation) matrix
correlation_matrix <- cor(norm.counts)
# View the correlation matrix
correlation_matrix[1:5, 1:5]

### Adjacancy matrix
# Calculate the adjacency matrix based on the chosen power (beta)
A = adjacency(norm.counts, power = soft_power)

##4B memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 2000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
cor <- temp_cor

#Module Eigengenes 
module_eigengenes <- bwnet$MEs
# get number of genes for each module
table(bwnet$colors)

#Plot dendogram and modules
plotDendroAndColors(bwnet$dendrograms[[1]], bwnet$colors,
                    "Modules",
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#module trait associations
# create traits file - binarize categorical variables
trait<-colData %>% 
  mutate(gender_bin = ifelse(grepl('female', gender), 1, 0)) %>% 
  select(gender_bin)

# binarize categorical variables
colData$type <- factor(colData$diagnosis, levels = c("Craniopharyngioma", "Medulloblastoma",   "Glioblastoma","Glioma", "Ependymoma","ATRT" ))

type.out <- binarizeCategoricalColumns(colData$type,
                                       includePairwise = FALSE,
                                       includeLevelVsAll = TRUE,
                                       minCount = 1)
colnames(type.out)

traits <- cbind(trait, type.out)
# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)
nSamples;nGenes

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr

module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
module.trait.corr.pvals

# visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% ### renames the columns for visisbility
  rename_all(~str_replace_all(., "data.", "")) %>%
  rename_all(~str_replace_all(., ".vs.all", "")) %>%
  rename_all(~str_replace_all(., "_bin", ""))


colnames <- names(heatmap.data)
# Rename the selected ones
colnames[1:8] <- c("Module1", "Module2", "Module3", "Module4", "Module5", "Module6", "Module7", "Module8")
# Assign the new names back to the dataframe
names(heatmap.data) <- colnames
cbbPalette <- c("#D55E00","#009E73")
par(cex.axis=0.8, cex.main=1.5)
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[11:15],
             y = names(heatmap.data)[1:8],
             cexLabX = .7,
             titleX = "tumor types",
             signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
             cexMain = 1,
             colFrame = "white",
             col = cbbPalette
)

#### Module mapping
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping%>%head() 

##### ADD Gene Names
#Intramodular analysis: Identifying driver genes
#Calculate the module membership and the associated p-values
#intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
module.membership.measure.pvals[1:9,1:10]

turquoise<-module.membership.measure.pvals["MEturquoise",]
hub_t<-as.data.frame(turquoise)%>%
  arrange(turquoise)%>%
  head()
ATRT_immune_res[rownames(hub_t),]

###Using the gene significance you can identify genes that have a high significance for trait of interest 
#Using the module membership measures you can identify genes with high module membership with module
#who has high correlation with Meduloblastoma.

# Calculate the gene significance and associated p-values
gene.signf.corr <- cor(norm.counts, traits$data.Medulloblastoma.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(5)

top2_MB<- Medulloblastoma_immune_res[c("ENSG00000033800.13", "ENSG00000115738.10"), ]
top2_MB

#Identifying top 5 module hubs
hub_genes <- function(module_color) {
  ME <- module.membership.measure.pvals[paste0("ME", module_color), ]
  module_genes <- names(bwnet$colors)[bwnet$colors == module_color]
  module_pvals <- ME[module_genes]
  hub_genes <- sort(module_pvals, decreasing = FALSE)
  top_5 <- hub_genes[1:5]
  top_5_df <- data.frame(ID = names(top_5), pvalue = top_5)
  top_5_df <- dplyr::inner_join(top_5_df, gene, by = "ID")
  return(top_5_df)
}
colors<-unique(bwnet$colors)[1:8]
for (i in 1:length(colors)) {
  print(bwnet$colors[i])
  print(hub_genes(colors[i]))
}

### END ################

