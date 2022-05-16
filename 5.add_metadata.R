# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(data.table)


# Define paths
path_to_data <- ("~/Documents/multiome_tonsil_Lucia/results/R_objects/6.tonsil_filtered_merged_all.rds")
path_to_multiome_metadata <- here::here("Documents/multiome_tonsil_Lucia/data/data/tonsil_atlas_metadata_multiome.csv")
path_to_scrublet_dfs <- here::here("Documents/multiome_tonsil_Lucia/2.doublet_detection/tmp/each_sample")
path_to_donor_metadata<-here::here("Documents/multiome_tonsil_Lucia/data/data/tonsil_atlas_donor_metadata.csv")
path_to_SCS_df <- here::here("Documents/multiome_tonsil_Lucia/results/tables/SCS_df.csv")
path_to_save <- ("~/Documents/multiome_tonsil_Lucia/results/R_objects/7.tonsil_filtered_merged_with_metadata.rds")


# Load data
tonsil_multiome_filtered_merged <- readRDS(path_to_data)

#become the first rowname as a column called "barcodes"
tonsil_multiome_filtered_merged@meta.data <- tibble::rownames_to_column(tonsil_multiome_filtered_merged@meta.data, "lib_name_barcode")
rn<-tonsil_multiome_filtered_merged@meta.data$lib_name_barcode
m1<-str_extract(rn, "BCLL\\_[\\d|\\d{2}]+\\_\\w\\_\\d")
tonsil_multiome_filtered_merged@meta.data$library_name<-m1

multiome_metadata <- read_csv(path_to_multiome_metadata)
donor_metadata <- read_csv(path_to_donor_metadata)
SCS_df <- read_csv(path_to_SCS_df)

scrublet_files <- list.files(path_to_scrublet_dfs)
sample_name_list<- scrublet_files
scrublet_files <- str_c(path_to_scrublet_dfs, scrublet_files, sep = "/")
all_scrublet <- purrr::map(str_c(scrublet_files,"data_frames/scrublet_doublet_prediction.csv", sep = "/"), read_csv)
names(all_scrublet)<-c(sample_name_list)

for (x in names(all_scrublet)){all_scrublet[[x]]$library_name<-x}
scrublet_df <- rbindlist(all_scrublet)
scrublet_df$lib_name_barcode<- paste(scrublet_df$library_name, scrublet_df$barcodes, sep="_")




# Include metadata
multiome_metadata_sub <- multiome_metadata %>%
  group_by(gem_id) %>%
  dplyr::filter(row_number(`gem_id`) == 1)
donor_ids <- multiome_metadata_sub$donor_id
names(donor_ids) <- multiome_metadata_sub$library_name
tonsil_multiome_filtered_merged@meta.data$donor_id <- donor_ids[tonsil_multiome_filtered_merged$library_name]
sex_vec <- donor_metadata$sex
age_vec <- donor_metadata$age
age_group_vec <- donor_metadata$age_group
hospital_vec <- donor_metadata$hospital
names(sex_vec) <- donor_metadata$donor_id
names(age_vec) <- donor_metadata$donor_id
names(age_group_vec) <- donor_metadata$donor_id
names(hospital_vec) <- donor_metadata$donor_id
tonsil_multiome_filtered_merged@meta.data$sex <- sex_vec[tonsil_multiome_filtered_merged$donor_id]
tonsil_multiome_filtered_merged@meta.data$age <- age_vec[tonsil_multiome_filtered_merged$donor_id]
tonsil_multiome_filtered_merged@meta.data$age_group <- age_group_vec[tonsil_multiome_filtered_merged$donor_id]
tonsil_multiome_filtered_merged@meta.data$hospital <- hospital_vec[tonsil_multiome_filtered_merged$donor_id]
tonsil_multiome_filtered_merged@meta.data$assay <- "multiome"
tonsil_multiome_filtered_merged@meta.data$S.Score<- SCS_df$S.Score
tonsil_multiome_filtered_merged@meta.data$G2M.Score<- SCS_df$G2M.Score
tonsil_multiome_filtered_merged@meta.data$Phase<- SCS_df$Phase


tonsil_multiome_filtered_merged@meta.data <- left_join(
  x = tonsil_multiome_filtered_merged@meta.data,
  y = scrublet_df,
  by = "lib_name_barcode")


##**change library_name column name** 
  
colnames(tonsil_multiome_filtered_merged@meta.data)[16]<-"library_name"
tonsil_multiome_filtered_merged@meta.data$library_name.y<-NULL

##**Make a rownames column** 
rownames(tonsil_multiome_filtered_merged@meta.data)<-tonsil_multiome_filtered_merged@meta.data$lib_name_barcode


write.csv(scrublet_df,here::here("Documents/multiome_tonsil_Lucia/results/scrublet_df.csv"))


# Save
saveRDS(tonsil_multiome_filtered_merged, path_to_save)




