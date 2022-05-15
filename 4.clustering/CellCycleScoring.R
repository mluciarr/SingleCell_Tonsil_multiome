library(Seurat)
library(dplyr)
library(Signac)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(plyr)
library(reshape2)
library(data.table)
library(GenomicRanges)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(hdf5r)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(magick)
library(knitr) 
library(kableExtra)
library(scater)
set.seed(123)
setwd("./Analysis")
library(devtools)

install_github("FloWuenne/scFunctions")

######Setup the Seurat obj

tonsil_multiome_QC <- readRDS("~/Documents/multiome_tonsil_Lucia/results/R_objects/1.tonsil_multiome_QC.rds")

normalize.list.SO <- function(seurat_object.list){
  
  seurat_object <- NormalizeData(seurat_object.list)
  
  return(seurat_object)
  
}

# CycleSocringCell

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

CycleScoringCell.list.SO <- function(seurat_object.list){
  
  seurat_object <- CellCycleScoring(seurat_object.list, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  return(seurat_object)
  
}


split_obj_function<- function(tonsil_filtered_RNA){
  for (i in seq_along(tonsil_filtered_RNA)){
    object<-tonsil_filtered_RNA[[i]]
    print(assign(paste("object",i,sep='_'), object, envir=.GlobalEnv))
  }
}  


remove_atac_assay <- function(seurat_object){
  DefaultAssay(seurat_object) <- "RNA"
  seurat_object@assays[["ATAC"]] <- NULL
    return(seurat_object)
}


# Normalizing the data

all_data<-lapply(tonsil_multiome_filtered, normalize.list.SO)

# CycleSocringCell

all_data_CSC<-lapply(all_data, CycleScoringCell.list.SO)



all_data_CSC_rna <- lapply(all_data_CSC, remove_atac_assay )




split_obj_function(all_data_CSC_rna)

merged<- merge(object_1, y = list(object_2,object_3,object_4,object_5,object_6,object_7,object_8,object_9,object_10,object_11), 
                          add.cell.ids = c(names(all_data_CSC) ))
rm(list=ls(pattern="object_"))

export_data_from_seurat(merged, output_dir = ".")

df<-as.data.frame(merged@meta.data)

write.csv(df,here::here("Documents/multiome_tonsil_Lucia/results/tables/SCS_df.csv"))

