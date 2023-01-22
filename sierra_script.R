sample_metadata <- read.table("sample_metadata.txt" , header = T , sep = ",") 
sample_metadata <- sample_metadata[sample_metadata$GEO_Accession..exp. != "GSM3516670",]

sample_metadata <- as.data.frame(cbind(sample_metadata$GEO_Accession..exp. , sample_metadata$Tissue , sample_metadata$STAGE))
colnames(sample_metadata) <- c("GEO" , "Tissue" , "Stage")
sample_metadata <- unique(sample_metadata)
sample_metadata$Source <- c("Tumor","Tumor","Met","Tumor","Normal","Tumor","Met",
                            "Tumor","Met","Tumor","Normal","Tumor","Normal","Normal","Met","Met")
geo_ids <- unique(sample_metadata$GEO)

reference.file <- "/media/huk/depo/references/GRCh38/GRCh38.108_ensembl.chr_filtered.gtf"

junction_files <- c()

for (i in 1:length(geo_ids)) {
  junction_files[i] <- paste("/media/huk/PirLab/gse123904/Cellranger_Out/",geo_ids[i],"/outs/possorted_genome_bam.bed" , sep = "")
}

bam_files <- c()

for (i in 1:length(geo_ids)) {
  bam_files[i] <- paste("/media/huk/PirLab/gse123904/Cellranger_Out/",geo_ids[i],"/outs/possorted_genome_bam.bam" , sep = "")
}

whitelist.bc.files <- c()

for (i in 1:length(geo_ids)) {
  whitelist.bc.files[i] <- paste("/media/huk/PirLab/gse123904/Cellranger_Out/",geo_ids[i],"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" , sep = "")
}


peak.output.files <- c()

for (i in 1:length(geo_ids)) {
  peak.output.files[i] <- paste("/home/huk/Documents/R/phd/GSE123904/sierra_peaks/",geo_ids[i],"_peaks.txt" , sep = "")
}

library(Sierra)


for (i in 1:length(geo_ids)) {
  FindPeaks(output.file = peak.output.files[i],   
            gtf.file = reference.file,           
            bamfile = bam_files[i],                
            junctions.file = junction_files[i],     
            ncores = 4)
}


peak.dataset.table = data.frame(Peak_file = peak.output.files,
                                Identifier = geo_ids, 
                                stringsAsFactors = FALSE)


peak.merge.output.file = "/home/huk/Documents/R/phd/GSE123904/sierra_peaks/GSE123904_merged_peaks.txt"

MergePeakCoordinates(peak.dataset.table, output.file = peak.merge.output.file, ncores = 5)


count.dirs <- c()


for (i in 1:length(geo_ids)) {
  count.dirs[i] <- paste("/home/huk/Documents/R/phd/GSE123904/",geo_ids[i],"_counts" , sep = "")
}


for (i in 1:length(geo_ids)) {
  CountPeaks(peak.sites.file = peak.merge.output.file, 
             gtf.file = reference.file,
             bamfile = bam_files[i], 
             whitelist.file = whitelist.bc.files[i],
             output.dir = count.dirs[i], 
             countUMI = TRUE, 
             ncores = 4)
}


out.dir <- "gse123904_peaks_aggregate"

AggregatePeakCounts(peak.sites.file = peak.merge.output.file,
                    count.dirs = count.dirs,
                    exp.labels = geo_ids,
                    output.dir = out.dir)


genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

AnnotatePeaksFromGTF(peak.sites.file = peak.merge.output.file, 
                     gtf.file = reference.file,
                     output.file = "gse123904_merged_peak_annotations.txt", 
                     genome = genome)


library(dplyr)
peak_annotations <- read.table("gse123904_merged_peak_annotations.txt", header = T, sep = "\t", 
                               row.names = 1, 
                               stringsAsFactors = FALSE) %>% filter(!is.na(gene_id))


peak.counts <- ReadPeakCounts(data.dir = out.dir)

