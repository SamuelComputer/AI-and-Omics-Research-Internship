# 1. Set Working Directory
getwd()
setwd("/home/samuel/Desktop/AI_Omics_Internship_2025/Module_I/")


  classify_gene <- function(logFC, padj) {
    if ((logFC > 1) & (padj < 0.05)) {
      return("Upregulated")
    } else if ((logFC < -1) & (padj < 0.05)) {
      return("Downregulated")
    } else {
      return("Not_Significant")
    }
  }

  
  input_dir <- "raw_data"
  output_dir <- "results"  
  
  # Create the folder if not exist
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    
  }
  
  file_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")
  
  # Prepare empty list to store results in R
  results_lists <- list()
  
  for (file_names in file_to_process) {
    cat("\nProcessing:", file_names, "\n")
  }
    
    input_file_path <- file.path(input_dir, file_names)
    
    # Import dataset
    data <- read.csv(input_file_path, header = TRUE)
    cat("File imported. Checking for missing values...\n")
    
    # Handling missing values
    if ("padj" %in% names(data)) {
      missing_count <- sum(is.na(data$padj))
      cat("Missing values in padj:", missing_count, "\n")
      
      # Replace missing padj with 1
      data$padj[is.na(data$padj)] <- 1
    }
    
    # Creating new column "status"
    data$status <- mapply(classify_gene, data$logFC, data$padj)
    
    #save results in R
    results_lists [[file_names]] <- data
    
    cat("Preview of data:\n")
    print(head(data, 10))
    View(data) 
    
    # Save processed file
    output_file_path <- file.path(output_dir, paste0("Processed_", file_names))
    write.csv(data, output_file_path, row.names = FALSE)
    cat("Processed file saved to:", output_file_path, "\n")
    
    
    # Summary counts
    summary_counts <- table(data$status)
    
    cat("\nSummary of gene classification:\n")
    print(summary_counts)
    