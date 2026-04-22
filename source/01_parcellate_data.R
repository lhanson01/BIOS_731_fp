library(here)
library(RNifti)
library(tidyverse)
library(ggnewscale)
library(cowplot)

#Load participants file #
participants <- read.csv(here("data/Meditation/participants.tsv"))

#### Load Data ####
data_path <- here("data", "Meditation")
parcellation_path <- here("data","Schaefer2018_300Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz")
#parcellation <- readNifti(parcellation_path)
sum(parcellation != 0)/length(parcellation)


ROIs <- unique(as.numeric((parcellation)))
mask_file <- file.path(data_path, "gm_mask_2mm.nii")
mask <- readNifti(mask_file)

anatomical_path <- file.path(data_path,"MNI152_T1_2mm.nii")
anatomical_data <- readNifti(anatomical_path)

plot_all_ROIs(
  template = anatomical_data,
  parcellation = parcellation,
  parc_title = "plot test",
  custom_slice = c(50,20,40)
)

#### Create vector of functional image file names ####
subject_paths <- file.path(data_path,
                           list.files(here("data", "Meditation"), pattern = "sub"))
n <- length(subject_images)


### generate n TxD matrices
subj_TxDs <- vector(mode = "list", length = 0)
for(i in 1:n){
  subj_TxDs[[paste("Subject", i)]] <-
    apply_parcellation(
      parcell_path = parcellation_path,
      funct_path = subject_paths[i]
    )
}

subj_ctrl <- 1:12
subj_trt <- 13:24

combined_ctrls <- bind_rows(subj_TxDs[1:12]) %>% filter(row_number() != 301)
combined_trts <- bind_rows(subj_TxDs[13:24]) %>% filter(row_number() != 301)

saveRDS(combined_ctrls,
        file = here::here("output", "combined_control_TxD.rds"))
saveRDS(combined_trts,
        file = here::here("output", "combined_trt_TxD.rds"))

con_mat_ctrl <- cor(t(combined_ctrls))
con_mat_trt <- cor(t(combined_trts))

plot_ctrl <- make_cor_heatmap(con_mat_ctrl, "Controls Connectivity")
plot_trt <- make_cor_heatmap(con_mat_trt, "Meditators Connectivity")

ggsave(file.path(here::here("output", "figures"), "control_heatmap.png"),
       plot_ctrl)

ggsave(file.path(here::here("output", "figures"), "med_heatmap.png"),
       plot_trt)



