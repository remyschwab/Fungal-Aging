library(DESeq2)
library(dplyr)
library(ggplot2)


read_col_data <- function(sample_sheet_path) {
  sample_sheet <- read.table(sample_sheet_path, header=TRUE)
  sample_sheet$age = factor(sample_sheet$age)
  sample_sheet$strain = factor(sample_sheet$strain)
  return(sample_sheet)
}


merge_count_files <- function(directory, pattern = "*.counts", count_column = 2, id_column = 1) {
  
  files <- list.files(path = directory, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop("No files found in the specified directory.")
  }
  
  count_list <- list()
  for (file in files) {
    data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    sample_name <- tools::file_path_sans_ext(basename(file))
    count_data <- data %>%
      select(!!id_column, !!count_column) %>%
      rename(TranscriptID = !!id_column, !!sample_name := !!count_column)
    count_list[[sample_name]] <- count_data
  }
  
  # Merge all data frames by the TranscriptID column
  merged_counts <- Reduce(function(x, y) full_join(x, y, by = "TranscriptID"), count_list)
  merged_counts[is.na(merged_counts)] <- 0
  
  # Change row names to values in the first column
  row.names(merged_counts) = merged_counts$TranscriptID
  merged_counts$TranscriptID = NULL
  
  return(merged_counts)
}

order_and_filter_results <- function(res, padj_threshold=0.05) {
  ordered = res[order(res$pvalue), ]
  filtered = ordered[which(ordered$padj<padj_threshold),]
  return(filtered)
}

plot_a_gene <- function(gene, dds) {
  plotCounts(dds, gene=gene, intgroup=c("strain", "age"), returnData=TRUE) %>%
    ggplot(aes(x=age, y=count, color=strain, group=strain)) +
    geom_point() + geom_line() +
    scale_y_log10() +
    theme_minimal()
}




## Summarize differential expression - Using the provided sample-level counts
## and sample sheet answer the following questions:
colData = read_col_data("calico_data_challenge_samplesheet.tsv")
countData = merge_count_files("RNAseq_counts/")
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~age + strain + strain:age)

## Which genes are differentially expressed based on age?
dds_age <- DESeq(dds, test = "LRT", reduced = ~ strain)
age_results = results(dds_age)

## Which genes are differentially expressed based on Genotype?
dds_genotype <- DESeq(dds, test = "LRT", reduced = ~ age)
genotype_results = results(dds_genotype)

## Which genes are differentially expressed based on genotype-specific aging?
dds_interact <- DESeq(dds, test = "LRT", reduced = ~ age + strain)
interact_results = results(dds_interact)

## Sanity Checks
## The counts of Knocked-out Genes should be lower in the corresponding strain
gene <- "YDR110W" # FOB1
plot_a_gene(gene, dds)
gene <- "YDL042C" # SIR2
plot_a_gene(gene, dds)
gene <- "YLR024C" # UBR2
plot_a_gene(gene, dds)
