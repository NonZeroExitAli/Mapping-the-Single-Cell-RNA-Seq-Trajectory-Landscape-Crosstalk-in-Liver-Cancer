library(tidyverse)
library(SingleCellExperiment)
library(celldex)
library(SingleR)
library(monocle3)
library(Seurat)
library(patchwork)

data_dir <- "E:/BioInfo_Project/FinaLProjectData/HCCT"
content_of_dir <- list.files(data_dir) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data_dir) #solving of error without gene.column 
scet1 = SingleCellExperiment(assays = list(counts = expression_matrix))

data_dir2 <- "E:/BioInfo_Project/FinaLProjectData/HCC02T"
content_of_dir2 <- list.files(data_dir2) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix2 <- Read10X(data_dir2) #solving of error without gene.column 
scet2 = SingleCellExperiment(assays = list(counts = expression_matrix2))

data_dir3 <- "E:/BioInfo_Project/FinaLProjectData/HCC03T"
content_of_dir3 <- list.files(data_dir3) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix3 <- Read10X(data_dir3) #solving of error without gene.column 
scet3 = SingleCellExperiment(assays = list(counts = expression_matrix3))


data_dir4 <- "E:/BioInfo_Project/FinaLProjectData/HCC04T"
content_of_dir4 <- list.files(data_dir4) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix4 <- Read10X(data_dir4) #solving of error without gene.column 
scet4 = SingleCellExperiment(assays = list(counts = expression_matrix4))



data_dir5 <- "E:/BioInfo_Project/FinaLProjectData/HCC05T"
content_of_dir5 <- list.files(data_dir5) # this line will show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix5 <- Read10X(data_dir5) #solving of error without gene.column 
scet5 = SingleCellExperiment(assays = list(counts = expression_matrix5))


scet <- cbind(scet1, scet2, scet3, scet4, scet5)



# Quality control (using mitochondrial genes)
library(scater)
mito <- grepl("^MT-", rownames(scet))
#grep <- return indices
# for grepl function <- function in scater package used to match the data using a string which needed and the object as inputs.
# return logical data (T, F)
summary(mito)
qcstats <- perCellQCMetrics(scet, subsets = list(Mito = mito))
#for subset argue <- A named list containing one or more vectors 
#(a character vector of feature names, a logical vector, or a numeric vector of indices),
#used to identify interesting feature subsets such as ERCC spike-in transcripts or mitochondrial genes.

filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
scet <- scet[, !filtered$discard]


# Normalization
library(scran)
#clusters <- quickCluster(scet)
#scet <- computeSumFactors(scet, clusters=clusters)
scet <- logNormCounts(scet)


# Feature selection
library(scran)
dec <- modelGeneVar(scet)
hvg <- getTopHVGs(dec, prop=0.1)


# PCA
library(scater)
set.seed(1234)
scet <- runPCA(scet, ncomponents=20, subset_row=hvg)
scet <- runTSNE(scet, dimred = 'PCA')


# Clustering
library(bluster)
colLabels(scet) <- clusterCells(scet, use.dimred='PCA',BLUSPARAM=NNGraphParam(cluster.fun="louvain"))   
#cllabels for work on colms(cells), object,indicates that clustering should be performed on the PCA scores
# specifies that the Louvain algorithm should be used for community detection 

# Visualization
scet <- runUMAP(scet, dimred = 'PCA') 
ug = plotUMAP(scet, color_by = 'label')
ug
tg=plotTSNE(scet, color_by = 'label')
tg
patchwork::wrap_plots(ug, tg, widths = 1, heights = 1, ncol = 1)

plotPCA(scet, color_by = 'label')

#markers detection (only if we will perform manual annotation)
markers <- findMarkers(scet, test.type="wilcox", direction="up", lfc=1)
marker.set <- markers[["20"]]
as.data.frame(marker.set[1:3,1:2])

marker.info.lfc <- scoreMarkers(scet, colLabels(scet), lfc=2)
chosen2 <- marker.info.lfc[["5"]] # another cluster for some variety.
chosen2 <- chosen2[order(chosen2$mean.AUC, decreasing=TRUE),]
chosen2[,c("self.average", "other.average", "mean.AUC")]
plotDots(scet, rownames(chosen2)[1:10], group="label")

plotExpression(scet, features=c("IGHA1"), x="label", colour_by="label")
plotDots(scet, rownames(chosen2)[1:10], group="label")

#traj
library(scater)
by.cluster <- aggregateAcrossCells(scet, ids=colLabels(scet))
centroids <- reducedDim(by.cluster, "PCA")
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst
line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")
dim(scet)
plotTSNE(scet, colour_by="label") + 
  geom_line(data=line.data, mapping=aes(x, y, group=edge))
