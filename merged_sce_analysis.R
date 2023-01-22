library(scater)
library(scran)
library(DropletUtils)
library(scuttle)
library(robustbase)
library(bluster)
library(Matrix)
library(dendextend)

load("~/Documents/R/phd/GSE123904/merged_sce.RData")
cell_info <- t(as.data.frame(str_split(merged_sce$cell_type , pattern = "_")))
merged_sce$cell_type <- cell_info[,1]
merged_sce$source <- cell_info[,2]
merged_sce$cell_label <- paste(cell_info[,1] , cell_info[,2] , sep = "_")

lib.lung <- librarySizeFactors(merged_sce)
summary(lib.lung)

hist(log10(lib.lung), xlab="Log10[Size factor]", col='grey80')

set.seed(125)

clust.lung <- quickCluster(merged_sce) 
table(clust.lung)

deconv <- calculateSumFactors(merged_sce, cluster=clust.lung)
summary(deconv)

plot(lib.lung, deconv, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16)
abline(a=0, b=1, col="red")

merged_sce <- computeSumFactors(merged_sce, cluster=clust.lung, min.mean=0.1)

merged_sce <- logNormCounts(merged_sce)
assayNames(merged_sce)

################### FEATURE SELECTION ##################

dec.lung <- modelGeneVar(merged_sce)

fit.lung <- metadata(dec.lung)

plot(fit.lung$mean, fit.lung$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.lung$trend(x), col="dodgerblue", add=TRUE, lwd=2)

dec.lung[order(dec.lung$bio, decreasing=TRUE),] 

hvg.lung.var <- getTopHVGs(dec.lung, prop = 0.1)
str(hvg.lung.var)

################### DIMENSION REDUCTION ###############

merged_sce <- denoisePCA(merged_sce , dec.lung ,
                                  subset.row = hvg.lung.var)

ncol(reducedDim(merged_sce, "PCA"))

plotReducedDim(merged_sce, dimred="PCA")

merged_sce <- runUMAP(merged_sce , dimred = "PCA")

plotReducedDim(merged_sce , dimred = "UMAP" , colour_by = "cell_type")
plotReducedDim(merged_sce , dimred = "UMAP" , colour_by = "source")

merged_sce <- runTSNE(merged_sce, dimred="PCA" , 
                               perplexity = 80 , num_threads = 8)

plotReducedDim(merged_sce , dimred = "TSNE" , colour_by = "cell_type")
plotReducedDim(merged_sce , dimred = "TSNE" , colour_by = "source")

######################## CLUSTERING ########################

lung.clusters <- clusterCells(merged_sce, use.dimred="PCA" ,
                              BLUSPARAM=SNNGraphParam(k=10, type="rank", cluster.fun="louvain"))

colLabels(merged_sce) <- lung.clusters

lung.clusters.info <- clusterCells(merged_sce, use.dimred="PCA" , full = TRUE ,
                                   BLUSPARAM=SNNGraphParam(k=10, type="rank", cluster.fun="louvain",
                                   ))


reducedDim(merged_sce, "force") <- igraph::layout_with_fr(lung.clusters.info$objects$graph)
plotReducedDim(merged_sce, colour_by="label", dimred="force")

plotReducedDim(merged_sce , dimred = "UMAP" , colour_by = "label" , text_by = "label")


