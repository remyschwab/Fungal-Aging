library(DESeq2)
library(dplyr)


## Summarize differential expression - Using the provided sample-level counts
## and sample sheet answer the following questions:

read_col_data <- function(sample_sheet_path) {
  sample_sheet <- read.table("WriteUp/calico_data_challenge_samplesheet.txt", header=TRUE)
  sample_sheet$age = factor(sample_sheet$age)
  sample_sheet$strain = factor(sample_sheet$strain)
  return(sample_sheet)
}


merge_count_files <- function(directory, pattern = "*.txt", count_column = 2, id_column = 1) {
  
  # List all files in the directory that match the pattern
  files <- list.files(path = directory, pattern = pattern, full.names = TRUE)
  
  # Check if there are files to process
  if (length(files) == 0) {
    stop("No files found in the specified directory.")
  }
  
  # Initialize an empty list to store data frames
  count_list <- list()
  
  # Iterate through each file and read the data
  for (file in files) {
    # Read the file into a data frame
    data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Extract sample name from the filename (can be customized)
    sample_name <- tools::file_path_sans_ext(basename(file))
    
    # Select the ID and count columns, renaming the count column to the sample name
    count_data <- data %>%
      select(!!id_column, !!count_column) %>%
      rename(TranscriptID = !!id_column, !!sample_name := !!count_column)
    
    # Add the data frame to the list
    count_list[[sample_name]] <- count_data
  }
  
  # Merge all data frames by the TranscriptID column
  merged_counts <- Reduce(function(x, y) full_join(x, y, by = "TranscriptID"), count_list)
  
  # Replace NAs with 0 (optional, depends on your use case)
  merged_counts[is.na(merged_counts)] <- 0
  
  # Change row names to values in the first column
  row.names(merged_counts) = merged_counts$TranscriptID
  merged_counts$TranscriptID = NULL
  
  # Return the merged count matrix
  return(merged_counts)
}


## Which genes are differentially expressed based on age?
dds_age <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~age)
dds_age <- DESeq(dds_age)
age_results <- results(dds_age, contrast=c("age", "40", "0"))

## Which genes are differentially expressed based on Genotype?
dds_strain <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~strain)
dds_strain <- DESeq(dds_strain)
strain_results_FOB1 <- results(dds_strain, contrast=c("strain", "DBY1200", "∆FOB1"))
strain_results_SIR2 <- results(dds_strain, contrast=c("strain", "DBY1200", "∆SIR2"))
strain_results_UBR2 <- results(dds_strain, contrast=c("strain", "DBY1200", "∆UBR2"))

## ## Which genes are differentially expressed based on genotype-specific aging?
dds_interact <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~age + strain + strain:age)
dds_interact <- DESeq(dds_interact)
