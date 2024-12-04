setwd('/lustre/home/mjf221/entropy_deconv')
library(tidyverse)
library(data.table)
library(ncdf4) #for opening netcdf files
### SYS ARGS ###

### OPEN BULK SEQUENCING DATA ###
K_matrix_config <- read_csv("K_matrices/K_matrix_config/cedric_0inflatedPoisResultsBased_Kmat_config_hg19_mode.csv")
betas <- list.files('cfDNA_files/loyfer_cfDNA', recursive=TRUE)
K_matrix_list <- list()
tissue_name_list <- list()
beta_matrix <- data.frame(nrows = 29152891) #if hg38
for (file in betas){
    fname <- paste0('cfDNA_files/loyfer_cfDNA/',file)
    N <- file.info(fname)$size
    df <- matrix(readBin(fname, "integer", N, size = 1, signed = FALSE), N / 2, 2, byrow=TRUE)
    
    mean_rd = round(mean(df[,2]), digits = -1)
    if (mean_rd > 100){
        mean_rd = 100
    }
    K_matrix_selected <- as.character(filter(K_matrix_config, mu_rd == mean_rd)[,3])
    K_matrix_list <- c(K_matrix_list, K_matrix_selected)

    beta_vec <- df[,1] / df[,2] 
    tissue_name <-  sub("^([^/]+)/.*$", "\\1", file)
    tissue_name_list <- c(tissue_name_list, tissue_name)
    beta_matrix <- cbind(beta_matrix, beta_vec)
}
bulk_seq_data <- as.data.frame(beta_matrix[,2:ncol(beta_matrix)])
colnames(bulk_seq_data) <- betas
rownames(bulk_seq_data) <- as.character(1:nrow(beta_matrix))


### CALCULATE CETYGO SCORE ###
predicted_cell_type_proportions <- read_csv("results/loyfer_realcfDNA_cedric_0inflatedPoisResultsBased_Kmat_config_hg19_mode_2024-06-13_1534-38.csv") #read in your deconvolution results file containing deconvoluted cell type proportions
cetygo_score_vec <- vector()
for (sample in seq(1:ncol(bulk_seq_data))){
    ### OPEN MODEL FILE (K MATRIX) ###

    K_matrix <- nc_open(paste0("K_matrices/", K_matrix_list[sample])) #cell type deconvolution model reference matrix: rows are cpg sites, columns are cell types, contents are predicted methylation betas
    reference_profile <- as.data.frame(ncvar_get(K_matrix, "__xarray_dataarray_variable__"))
    colnames(reference_profile) <- predicted_cell_type_proportions$ct
    rownames(reference_profile) <- ncvar_get(K_matrix, "gen")

    ### CALCULATING PREDICTED BULK METHYLATION PROFILE ###
    predicted_cell_type_proportions_t <- t(predicted_cell_type_proportions[,2:ncol(predicted_cell_type_proportions)])
    colnames(predicted_cell_type_proportions_t) <- predicted_cell_type_proportions$ct 

    product <- data.frame(nrow=nrow(reference_profile))
    for(column in seq(ncol(predicted_cell_type_proportions_t))){
        cell_type_k_product <- predicted_cell_type_proportions_t[sample,column] * reference_profile[,column]
        product <- cbind(product, cell_type_k_product)
    }
    predicted_bulk_profile <- rowSums(product[2:ncol(product)])

    ### CALCULATE CETYGO SCORE ###
    observed_bulk_profile <- beta_matrix[,2:ncol(beta_matrix)]
    observed_bulk_profile <- observed_bulk_profile[,sample]
    observed_bulk_profile_filtered <- observed_bulk_profile[as.numeric(rownames(reference_profile))]

    intermediate <- (observed_bulk_profile_filtered - predicted_bulk_profile)^2
    intermediate_noNaN <- intermediate[!is.na(intermediate)]
    cetygo_score <- sqrt(mean(intermediate_noNaN)) 
    cetygo_score_vec <- c(cetygo_score_vec, cetygo_score)
}

cetygo_results <- data.frame(sample = colnames(bulk_seq_data), CETYGO_score = cetygo_score_vec)
write.table(cetygo_results, "bactydeMode_loyfer_real_cfDNA_cetygo_results.csv", sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)