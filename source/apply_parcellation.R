apply_parcellation <- function(parcell_path, labels_path, funct_path){
  print("New Image!")
  funct_data <- RNifti::readNifti(funct_path)
  parcell <- RNifti::readNifti(parcell_path)
  ROI_table <- read.csv(labels_path)
  ROIs <- ROI_table$node_index
  names(ROIs) <- ROI_table$node_label
  n_ROI <- length(ROIs)
  nt <- dim(funct_data)[4]
  v <- prod(dim(funct_data)[1:3])
  func_vxt <- matrix(funct_data, nrow = v)
  ROI_bold_ts_mat <- matrix(data = NA, nrow = n_ROI, ncol = nt)
  for (ROI in ROIs) {
    ROI_vector_indeces <- which(parcell == ROI)
    ROI_bold_ts_mat[ROI,] <- apply(func_vxt[ROI_vector_indeces,], 2, mean)
  }
  return(ROI_bold_ts_mat)
}

# slow_parcellation_flipped <- function(parcell, funct_data){
#   ROIs <- unique(as.numeric((parcell)))
#   n_ROI <- length(ROIs)
#   nt <- dim(funct_data)[4]
#   ROI_bold_ts_mat <- matrix(data = NA, nrow = n_ROI, ncol = nt)
#   for (t in 1:nt) {
#     for (ROI in ROIs){
#       ROI_vector_indeces <- which(parcell == ROI)
#       ROI_bold_ts_mat[ROI, t] <- mean(funct_data[,,,t][ROI_vector_indeces])
#     }
#   }
#   return(ROI_bold_ts_mat)
# }
