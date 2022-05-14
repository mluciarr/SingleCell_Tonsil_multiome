library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggpubr)

set.seed(173)


filtered.data<-readRDS("~/Documents/multiome_tonsil_Lucia/results/R_objects/2.tonsil_multiome_filtered.rds")


peaks_quantification <- function(seurat_filtered, peaks = combined.peaks){ 
  counts <- FeatureMatrix(
    fragments = Fragments(seurat_filtered),
    features = combined.peaks,
    cells = colnames(seurat_filtered))
  seurat_filtered[["peaks"]] <- CreateChromatinAssay(
    counts, 
    genome = "hg38",
    fragments = Fragments(seurat_filtered),
    annotation = annotation)
  
  return(seurat_filtered)
}

default_assays <- function(seurat_object){
  DefaultAssay(seurat_object) <- "ATAC"
  return(seurat_object)
}

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,standard.chromosomes = T)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# Create a unified set of peaks
tonsil_data_filtered <- lapply(filtered.data, default_assays)
combined.peaks <- UnifyPeaks(object.list = tonsil_data_filtered, mode = "reduce") 

# Keeping only standard chromosomes
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = 'coarse')
peakwidths <- width(combined.peaks)

p<-ggviolin(peakwidths,add = "boxplot",fill = "gray") + scale_y_log10() + 
  geom_hline(yintercept = c(20,10000), linetype='dashed', col = 'black')

combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
combined.peaks

tonsil_filtered_unified_peaks <- lapply(tonsil_data_filtered, peaks_quantification)

# Removing ATAC assay to reduce the size of the final seurat object,

remove_assays <- function(seurat_object){
  DefaultAssay(seurat_object) <- "peaks"
  seurat_object@assays[["ATAC"]] <- NULL
  
  RenameAssays(object = seurat_object, peaks = 'ATAC')
  return(seurat_object)
}

tonsil_filtered_unified_peaks <- lapply(tonsil_filtered_unified_peaks, remove_assays)
saveRDS(tonsil_filtered_unified_peaks,"../results/R_objects/4.tonsil_filtered_unified_peaks.rds")


#Mergin suerat objects
object1 <- tonsil_filtered_unified_peaks[[1]]
object2 <- tonsil_filtered_unified_peaks[[2]]
object3 <- tonsil_filtered_unified_peaks[[3]]
object4 <- tonsil_filtered_unified_peaks[[4]]
object5 <- tonsil_filtered_unified_peaks[[5]]
object6 <- tonsil_filtered_unified_peaks[[6]]
object7 <- tonsil_filtered_unified_peaks[[7]]
object8 <- tonsil_filtered_unified_peaks[[8]]
object9 <- tonsil_filtered_unified_peaks[[9]]
object10 <- tonsil_filtered_unified_peaks[[10]]
object11 <- tonsil_filtered_unified_peaks[[11]]

merged <- merge(object1, y = list(object2, object3, object4, object5, object6,
                                  object7, object8,object9,object10,object11), 
                add.cell.ids = c("BCLL_15_T_1", "BCLL_2_T_3", "BCLL_9_T_1",  "BCLL_2_T_1",
                                 "BCLL_14_T_1", "BCLL_9_T_2",  "BCLL_2_T_2",  "BCLL_15_T_2", "BCLL_8_T_1", 
                                 "BCLL_14_T_2", "BCLL_8_T_2"))
saveRDS(tonsil_filtered_unified_peaks,"../results/R_objects/4.tonsil_filtered_unified_peaks_merged.rds")

# Create matrices
rna_mat_bcll_2 <- filtered.data[["RNA"]]@counts[, bcll_2_cells]
rna_mat_non_bcll_2 <- filtered.data[["RNA"]]@counts[, non_bcll_2_cells]
atac_mat_bcll_2 <- filtered.data[["peaks"]]@counts[, bcll_2_cells]
atac_mat_non_bcll_2 <- filtered.data[["peaks"]]@counts[, non_bcll_2_cells]



# Save
dir.create(path_to_tmp_dir, showWarnings = FALSE)
mats_list <- list(
  rna_sparse_matrix_with_BCLL_2 = rna_mat_bcll_2,
  rna_sparse_matrix_without_BCLL_2 = rna_mat_non_bcll_2,
  atac_sparse_matrix_with_BCLL_2 = atac_mat_bcll_2,
  atac_sparse_matrix_without_BCLL_2 = atac_mat_non_bcll_2
)
for (x in names(mats_list)) {
  path_to_subdir <- str_c(path_to_tmp_dir, x, sep = "")
  dir.create(path_to_subdir, showWarnings = FALSE)
  path_save_mat <- str_c(path_to_subdir, "matrix.mtx", sep = "/")
  path_save_features <- str_c(path_to_subdir, "features.tsv", sep = "/")
  path_save_cell_barcodes <- str_c(path_to_subdir, "barcodes.tsv", sep = "/")
  writeMM(mats_list[[x]], path_save_mat)
  write(x = rownames(mats_list[[x]]), file = path_save_features)
  write(x = colnames(mats_list[[x]]), file = path_save_cell_barcodes)
}
