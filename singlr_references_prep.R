library(SingleR)

##### SINGLR REFERENCE #####

ref_info <- read.table("Cell_Info.tsv" , sep = "\t" , header = T)
ref_data <- read.table("DF_ALL.tsv" , sep = "\t" , header = T , row.names = 1)
load("~/Documents/R/phd/GSE123904/data_genes.RData")
common_genes_between_data_and_ref <- intersect(common_genes , rownames(ref_data))

tumor_ref <- ref_data[,ref_info$Meta.Source == "TUMOR"]
tumor_ref <- tumor_ref[common_genes_between_data_and_ref,]
colnames(tumor_ref) <- gsub("\\..*" , "" , colnames(tumor_ref))
tumor_ref_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = tumor_ref))
tumor_ref_sce <- scuttle::logNormCounts(tumor_ref_sce)
tumor_ref <- SingleCellExperiment::logcounts(tumor_ref_sce)

met_ref <- ref_data[,ref_info$Meta.Source == "MET"]
met_ref <- met_ref[common_genes_between_data_and_ref,]
colnames(met_ref) <- gsub("\\..*" , "" , colnames(met_ref))
met_ref_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = met_ref))
met_ref_sce <- scuttle::logNormCounts(met_ref_sce)
met_ref <- SingleCellExperiment::logcounts(met_ref_sce)

normal_ref <- ref_data[,ref_info$Meta.Source == "NOR"]
normal_ref <- normal_ref[common_genes_between_data_and_ref,]
colnames(normal_ref) <- gsub("\\..*" , "" , colnames(normal_ref))
normal_ref_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = normal_ref))
normal_ref_sce <- scuttle::logNormCounts(normal_ref_sce)
normal_ref <- SingleCellExperiment::logcounts(normal_ref_sce)

rm(ref_data, tumor_ref_sce , met_ref_sce , normal_ref_sce)


trained_met_ref <- trainSingleR(met_ref, colnames(met_ref))
trained_normal_ref <- trainSingleR(normal_ref, colnames(normal_ref))
trained_tumor_ref <- trainSingleR(tumor_ref, colnames(tumor_ref))

save(trained_met_ref , trained_normal_ref , trained_tumor_ref, common_genes_between_data_and_ref , file = "singler_references.RData")