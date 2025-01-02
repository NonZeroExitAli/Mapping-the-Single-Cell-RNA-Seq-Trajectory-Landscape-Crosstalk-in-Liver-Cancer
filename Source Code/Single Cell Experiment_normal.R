library(tidyverse)
library(SingleCellExperiment)
library(celldex)
library(SingleR)
library(Seurat)
library(patchwork)


data_dir <- "E:/BioInfo_Project/FinaLProjectData/HCC01P"
content_of_dir <- list.files(data_dir) # th is line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data_dir) #solving of error without gene.column 
scen1 = SingleCellExperiment(assays = list(counts = expression_matrix))


data_dir2 <- "E:/BioInfo_Project/FinaLProjectData/HCC02P"
content_of_dir2 <- list.files(data_dir2) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix2 <- Read10X(data_dir2) #solving of error without gene.column 
scen2 = SingleCellExperiment(assays = list(counts = expression_matrix2))


data_dir3 <- "E:/BioInfo_Project/FinaLProjectData/HCC03P"
content_of_dir3 <- list.files(data_dir3) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix3 <- Read10X(data_dir3) #solving of error without gene.column 
scen3 = SingleCellExperiment(assays = list(counts = expression_matrix3))


data_dir4 <- "E:/BioInfo_Project/FinaLProjectData/HCC04P"
content_of_dir4 <- list.files(data_dir4) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix4 <- Read10X(data_dir4) #solving of error without gene.column 
scen4 = SingleCellExperiment(assays = list(counts = expression_matrix4))



data_dir5 <- "E:/BioInfo_Project/FinaLProjectData/HCC05P"
content_of_dir5 <- list.files(data_dir5) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix5 <- Read10X(data_dir5) #solving of error without gene.column 
scen5 = SingleCellExperiment(assays = list(counts = expression_matrix5))


scenn <- cbind(scen1, scen2, scen3, scen4, scen5)



# Quality control (using mitochondrial genes)
library(scater)
mito <- grepl("^MT-", rownames(scenn))
#grep <- return indices
# for grepl function <- function in scater package used to match the data using a string which needed and the object as inputs.
# return logical data (T, F)
summary(mito)
qcstats <- perCellQCMetrics(scenn, subsets = list(Mito = mito))
#for subset argue <- A named list containing one or more vectors 
#(a character vector of feature names, a logical vector, or a numeric vector of indices),
#used to identify interesting feature subsets such as ERCC spike-in transcripts or mitochondrial genes.

filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
scenn <- scenn[, !filtered$discard]


# Normalization
library(scran)
#clusters <- quickCluster(scenn)
#scenn <- computeSumFactors(scenn, clusters=clusters)
scenn <- logNormCounts(scenn)


# Feature selection
library(scran)
dec <- modelGeneVar(scenn)
hvg <- getTopHVGs(dec, prop=0.1)


# PCA
library(scater)
set.seed(1234)
scenn <- runPCA(scenn, ncomponents=20, subset_row=hvg)
scenn <- runTSNE(scenn, dimred = 'PCA')


# Clustering
library(bluster)
colLabels(scenn) <- clusterCells(scenn, use.dimred='PCA',BLUSPARAM=NNGraphParam(cluster.fun="louvain"))   
#cllabels for work on colms(cells), object,indicates that clustering should be performed on the PCA scores
# specifies that the Louvain algorithm should be used for community detection 


# Visualization
scenn <- runUMAP(scenn, dimred = 'PCA') 
plotUMAP(scenn, color_by = 'label')
plotTSNE(scenn, color_by = 'label')

#markers detection (only if we will perform manual annotation)
markers <- findMarkers(scenn, test.type="wilcox", direction="up", lfc=1)
marker.set <- markers[["1"]]
as.data.frame(marker.set[1:3,1:2])


#Annotation (automated annotation)
library(celldex)
ref <- HumanPrimaryCellAtlasData()

library(SingleR)
pred <- SingleR(test=scenn, ref=ref, labels=ref$label.main)
table(pred$labels)

labels <- pred$labels
colData(scenn)$labels <- labels


# Visualization of annotatid cells

#automated annotation visualization
ug = plotUMAP(scenn, color_by = 'labels', text_by = "labels")
tg = plotTSNE(scenn, color_by = 'labels', text_by = "labels")
tg


markers <- findMarkers(scenn, test.type="wilcox", direction="up", lfc=1)
marker.set <- markers[["20"]]
as.data.frame(marker.set[1:3,1:2])

marker.info.lfc <- scoreMarkers(scenn, colLabels(scenn), lfc=2)
chosen2 <- marker.info.lfc[["5"]] # another cluster for some variety.
chosen2 <- chosen2[order(chosen2$mean.AUC, decreasing=TRUE),]
chosen2[,c("self.average", "other.average", "mean.AUC")]
plotDots(scenn, rownames(chosen2)[1:10], group="label")

plotExpression(scenn, features=c("IGHA1"), x="label", colour_by="label")
plotDots(scenn, rownames(chosen2)[1:10], group="label")


#traj
library(scater)
by.cluster <- aggregateAcrossCells(scenn, ids=colLabels(scenn))
centroids <- reducedDim(by.cluster, "PCA")
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst
line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")
dim(scenn)
plotTSNE(scenn, colour_by="label") + 
  geom_line(data=line.data, mapping=aes(x, y, group=edge))

