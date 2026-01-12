gdc <- readLines("gdc_472_HNSC.txt")
gdc <- gdc[-c(473:474)]
gdc_sample <- matrix(gdc, nrow = length(gdc), ncol =2)
patient_id <- paste0("T", c(1:472))
gdc_sample[,2] <- patient_id 
colnames(gdc_sample) <- c("file_id", "patient_id")
file_paths <- list.files(path = ".", pattern = "\\.tsv$", full.name = TRUE, recursive = TRUE)
tsv_list <- lapply(file_paths, readLines)
cleaner <- function(x) {x <- x[-c(1,3,4,5,6)]}
tumor_data <- lapply(tsv_list, cleaner)
convert_to_dataframe <- function(lines) {                                                                             
     header <- unlist(strsplit(lines[1], "\t"))                                                                        
     data_lines <- lines[-1]                                                                                           
        data <- do.call(rbind, lapply(data_lines, function(x) {                                                        
               unlist(strsplit(x, "\t"))                                                                               
         }))                                                                                                           
   df <- as.data.frame(data, stringsAsFactors = FALSE)                                                                 
   colnames(df) <- header                                                                                              
   return(df)                                                                                                          
}
tumor_df <- lapply(tumor_data, convert_to_dataframe)
selector <- function(x) {x <- x[which(x$gene_type == "protein_coding"),]}
tumor_df_pc <- lapply(tumor_df, selector)
selector2 <- function(x) {x <- x[,c(2,7)]}
tumor_df_pc_tpm <- lapply(tumor_df_pc, selector2)
library(dplyr)
distinct_tumor <- lapply(tumor_df_pc_tpm, function(df) distinct(df, gene_name, .keep_all = TRUE))
tumor <- Reduce(function(x, y) inner_join(x, y, by = "gene_name"), distinct_tumor)
colnames(tumor) <- c("gene_name", patient_id)
write.csv(tumor, file = "tumor_tpm.csv")
write.csv(gdc_sample, file = "gdc_file_id_to_patient_id.csv")