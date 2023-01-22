library(scater)
library(scran)
library(DropletUtils)
library(scuttle)
library(robustbase)
library(bluster)
library(Matrix)
library(dendextend)
library(devtools)

sample_metadata <- read.table("sample_metadata.txt" , header = T , sep = ",") 
sample_metadata <- sample_metadata[sample_metadata$GEO_Accession..exp. != "GSM3516670",]

sample_metadata <- as.data.frame(cbind(sample_metadata$GEO_Accession..exp. , sample_metadata$Tissue , sample_metadata$STAGE))
colnames(sample_metadata) <- c("GEO" , "Tissue" , "Stage")
sample_metadata <- unique(sample_metadata)
sample_metadata$Source <- c("Tumor","Tumor","Met","Tumor","Normal","Tumor","Met",
                            "Tumor","Met","Tumor","Normal","Tumor","Normal","Normal","Met","Met")

geo_ids <- unique(sample_metadata$GEO)
data_dirs <- list()

for (i in 1:length(geo_ids)) {
  data_dirs[i] <- paste("/media/huk/PirLab/gse123904/Cellranger_Out/",geo_ids[i],"/outs/filtered_feature_bc_matrix" , sep = "")
}

datasets <- list()

for (i in 1:length(geo_ids)) {
  datasets[[i]] <- read10xCounts(data_dirs[i] , col.names = TRUE)
}

for (i in 1:length(datasets)){
  rownames(datasets[[i]]) <- uniquifyFeatureNames(rowData(datasets[[i]])$ID,rowData(datasets[[i]])$Symbol)
}

for (i in 1:length(datasets)){
  is.mito <- grepl("^MT-" , rownames(datasets[[i]]))
  qc_df <- perCellQCMetrics(datasets[[i]] , subsets=list(Mito=is.mito))
  filters <- perCellQCFilters(qc_df,sub.fields=c("subsets_Mito_percent"))
  colData(datasets[[i]]) <- cbind(colData(datasets[[i]]) , qc_df)
  datasets[[i]]$discard <- filters$discard
  datasets[[i]] <- addPerCellQCMetrics(datasets[[i]], subsets=list(Mito=is.mito))
}



for (i in 1:length(datasets)){
  png(file = paste("qual_plots/",geo_ids[i],"_quality_plots.png",sep="") ,width = 1280 ,height = 720)
  gridExtra::grid.arrange(
    plotColData(datasets[[i]], y="sum", colour_by="discard") + 
      scale_y_log10() + ggtitle("Total count") + xlab(geo_ids[i]),
    plotColData(datasets[[i]], y="detected", colour_by="discard" ) +
      scale_y_log10() + ggtitle("Detected features") + xlab(geo_ids[i]),
    plotColData(datasets[[i]], y="subsets_Mito_percent", 
              colour_by="discard") +  ggtitle("Mito percent") + xlab(geo_ids[i]),
    ncol=1
)
  dev.off()
}



for (i in 1:length(datasets)){
  datasets[[i]] <- datasets[[i]][,!datasets[[i]]$discard]
  lib.lung <- librarySizeFactors(datasets[[i]])
  clust.lung <- quickCluster(datasets[[i]]) 
  deconv <- calculateSumFactors(datasets[[i]], cluster=clust.lung)
  datasets[[i]] <- computeSumFactors(datasets[[i]], cluster=clust.lung, min.mean=0.1)
  datasets[[i]] <- logNormCounts(datasets[[i]])
}
  

rm(filters , qc_df , is.mito , clust.lung , lib.lung , deconv)

common_genes <- intersect(rownames(datasets[[1]]) , rownames(datasets[[2]]))

for (i in 3:length(datasets)) {
  common_genes <- intersect(common_genes , rownames(datasets[[i]]))
}

for (i in 1:length(datasets)) {
  datasets[[i]] <- datasets[[i]][common_genes,]
}


############################# CELL META DATA ##########################3

sample_ids <- c(rep(geo_ids[1] , ncol(datasets[[1]])) , rep(geo_ids[2] , ncol(datasets[[2]])))

for (i in 3:length(datasets)) {
  sample_ids <- c(sample_ids , rep(geo_ids[i] , ncol(datasets[[i]])))
}

sample_tissue <- c(rep(sample_metadata$Tissue[1] , ncol(datasets[[1]])) , rep(sample_metadata$Tissue[2] , ncol(datasets[[2]])))

for (i in 3:length(datasets)) {
  sample_tissue <- c(sample_tissue , rep(sample_metadata$Tissue[i] , ncol(datasets[[i]])))
}

sample_stage <- c(rep(sample_metadata$Stage[1] , ncol(datasets[[1]])) , rep(sample_metadata$Stage[2] , ncol(datasets[[2]])))

for (i in 3:length(datasets)) {
  sample_stage <- c(sample_stage , rep(sample_metadata$Stage[i] , ncol(datasets[[i]])))
}


##############################################################################
library(SingleR)
load("~/Documents/R/phd/GSE123904/singler_references.RData")

dataset_labels <- list()

for (i in 1:length(datasets)) {
  if (sample_metadata$Source[i] == "Tumor") {
    dataset_labels[[i]] <- classifySingleR(test = datasets[[i]],trained = trained_tumor_ref)$labels
  }
  else if (sample_metadata$Source[i] == "Met") {
    dataset_labels[[i]] <- classifySingleR(test = datasets[[i]],trained = trained_met_ref)$labels
  }
  else if (sample_metadata$Source[i] == "Normal") {
    dataset_labels[[i]] <- classifySingleR(test = datasets[[i]],trained = trained_normal_ref)$labels
  }
}

merged_labels <- c(paste(dataset_labels[[1]],sample_metadata$Source[1],sep = "_") , 
                   paste(dataset_labels[[2]],sample_metadata$Source[2],sep = "_"))

for (i in 3:length(datasets)) {
  merged_labels <- c(merged_labels , paste(dataset_labels[[i]],sample_metadata$Source[i],sep = "_"))
}

################### MERGING COUNTS ######################

merged_counts <- cbind(counts(datasets[[1]]) , counts(datasets[[2]]))

for (i in 3:length(datasets)) {
  merged_counts <- cbind(merged_counts , counts(datasets[[i]]))
}

merged_sce <- SingleCellExperiment(assays = list(counts = merged_counts))

merged_sce$sample_ids <- sample_ids
merged_sce$sample_tissue <- sample_tissue
merged_sce$sample_stage <- sample_stage
merged_sce$cell_type <- merged_labels

save(merged_sce , file = "merged_sce.RData")

