
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(Repitools)

set.seed(123)

# Peak calling 

path_to_obj <- ("~/Documents/multiome_tonsil_Lucia/results/R_objects/14.tonsil_multiome_bcells_without_cluster4_doublets_normalized.rds")
path_to_save<- ("~/Documents/multiome_tonsil_Lucia/results/R_objects/14.1.tonsil_multiome_bcells_without_cluster4_doublets_normalized_Linkpeaks.rds")
path_to_save_normalized<- ("~/Documents/multiome_tonsil_Lucia/results/R_objects/14.2.tonsil_multiome_bcells_without_cluster4_doublets_normalized_Linkpeaks_normalized.rds")



## Load data

tonsil_wnn_bcell <- readRDS(path_to_obj)
#tonsil_wnn_bcell_peakcall<-readRDS(path_to_save)

DefaultAssay(tonsil_wnn_bcell) <- "peaks"


# call peaks using MACS2
peaks <- CallPeaks(tonsil_wnn_bcell, macs2.path = "~/Documents/python3/bin/macs2",group.by="wsnn_res.0.075")
write_csv(df, "~/Documents/multiome_tonsil_Lucia/results/tables/13.df_MACS_annotation_level_1.csv")
saveRDS(peaks, "~/Documents/multiome_tonsil_Lucia/results/tables/13.MACS_annotation_level_1.rds")
df <- annoGR2DF(peaks)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
saveRDS(peaks, "~/Documents/multiome_tonsil_Lucia/results/tables/13.peaks_annotation_level_1_subset.rds")

# quantify counts in each peak:
## Construct a feature x cell matrix from a genomic fragments file


macs2_counts <- FeatureMatrix(
  fragments = Fragments(tonsil_wnn_bcell),
  features = peaks,
  cells = colnames(tonsil_wnn_bcell)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

fragments<-Fragments(tonsil_wnn_bcell[["peaks"]])

tonsil_wnn_bcell[["ATAC"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragments,
  annotation = annotation
)


#Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
tonsil_wnn_bcell_peakcall<-tonsil_wnn_bcell
DefaultAssay(tonsil_wnn_bcell_peakcall) <- "ATAC"
tonsil_wnn_bcell_peakcall <- RunTFIDF(tonsil_wnn_bcell_peakcall)
tonsil_wnn_bcell_peakcall <- FindTopFeatures(tonsil_wnn_bcell_peakcall, min.cutoff = "q0")
tonsil_wnn_bcell_peakcall <- RunSVD(tonsil_wnn_bcell_peakcall)


saveRDS(tonsil_wnn_bcell, path_to_save)
saveRDS(tonsil_wnn_bcell_peakcall, path_to_save_normalized)
