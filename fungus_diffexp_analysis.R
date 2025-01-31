library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)


read_col_data <- function(sample_sheet_path) {
  sample_sheet <- read.table(sample_sheet_path, header=TRUE)
  sample_sheet$age = factor(sample_sheet$age)
  sample_sheet$strain = factor(sample_sheet$strain)
  row.names(sample_sheet) = sample_sheet$name
  return(sample_sheet)
}

merge_count_files <- function(directory, sample_sheet, count_column = 2, id_column = 1) {
  count_list <- list()
  for(i in 1:nrow(sample_sheet)) {
    count_file = sample_sheet$filename[i]
    data <- read.table(paste(directory, count_file, sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    sample_name <- sample_sheet$name[i]
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

plot_gene_expression <- function(dataset, standard_name, systematic_name) {
  df <- plotCounts(dds, gene=systematic_name, intgroup="strain", returnData=TRUE)
  ggplot(df, aes(x=strain, y=count, color=strain)) +
    geom_jitter(width=0.2, size=2) +
    geom_boxplot(alpha=0.5) +
    scale_y_log10() +  # Log-scale for better visualization
    labs(title=paste("Expression of", standard_name, "across strains"), y="Normalized Counts") +
    theme_minimal()
}



## Summarize differential expression
## Using the provided sample-level counts and sample sheet answer the following questions:
colData = read_col_data("calico_data_challenge_samplesheet.tsv")
countData = merge_count_files("RNAseq_counts/", colData)
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~age + strain + strain:age)

## Sanity Checks
## The counts of Knocked-out Genes should be lower in the corresponding strain
plot_gene_expression(dds, "FOB1", "YDR110W")
plot_gene_expression(dds, "SIR2", "YDL042C")
plot_gene_expression(dds, "UBR2", "YLR024C")

## Check for Sample mislabels
all(rownames(colData) == colnames(countData))
# colData <- colData[match(colnames(countData), rownames(colData)), ]
all(colData(dds)$name == colnames(counts(dds)))

vsd <- vst(dds)
pheatmap(cor(assay(vsd)), annotation_col = colData[,c("age","strain","batch")])
plotPCA(vsd, intgroup="age")
plotPCA(vsd, intgroup="strain")
plotPCA(vsd, intgroup="batch")

## LRT Version RIDICULOUS ADJUSTED PVALUES
## Which genes are differentially expressed based on age?
dds_age <- DESeq(dds, test = "LRT", reduced = ~ strain)
results_age = results(dds_age, alpha = 0.05)

## Which genes are differentially expressed based on Genotype?
dds_genotype <- DESeq(dds, test = "LRT", reduced = ~ age)
results_genotype = results(dds_genotype, alpha = 0.05)

## Which genes are differentially expressed based on genotype-specific aging?
dds_interact <- DESeq(dds, test = "LRT", reduced = ~ age + strain)
results_interact = results(dds_interact, alpha = 0.05)

## Wald Test Version
# Age
wald_age <- DESeq(dds)
results_wald_age <- results(wald_age, contrast = c("age", "40", "0"))
# Genotype
wald_genotype

