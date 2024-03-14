require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("mclust")
source("R/Analysis.R") # load after cyclemix
require("stringr")
library("Seurat")
library("biomaRt")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
options(max.print = 20)

convert_to_ensembl <- function(gene_symbols) {
    # Convert gene symbols to Ensembl IDs
    ensembl_ids_df <- select(org.Hs.eg.db, keys = gene_symbols, columns = "ENSEMBL", keytype = "SYMBOL")

    # Remove rows with NA ENSEMBL
    ensembl_ids_df <- ensembl_ids_df[!is.na(ensembl_ids_df$ENSEMBL), ]

    # If there are multiple ENSEMBLs for a single SYMBOL, keep the first one
    ensembl_ids_df <- ensembl_ids_df[!duplicated(ensembl_ids_df$SYMBOL), ]

    # # Create a named vector of ENSEMBLs with SYMBOLs as names
    # ensembl_ids <- setNames(ensembl_ids_df$ENSEMBL, ensembl_ids_df$SYMBOL)

    # # Replace the gene symbols with Ensembl IDs
    # ensembl_ids <- ensembl_ids[!is.na(gene_symbols)]

    return(ensembl_ids_df$ENSEMBL)
}

convert_to_symbols <- function(ensembl_ids) {
    # Convert Ensembl IDs to gene symbols
    gene_symbols_df <- select(org.Hs.eg.db, keys = ensembl_ids, columns = "SYMBOL", keytype = "ENSEMBL")

    # If there are multiple SYMBOLs for a single ENSEMBL, keep the first one
    gene_symbols_df <- gene_symbols_df[!duplicated(gene_symbols_df$ENSEMBL), ]

    # Create a named vector of SYMBOLs with ENSEMBLs as names
    gene_symbols <- setNames(gene_symbols_df$SYMBOL, gene_symbols_df$ENSEMBL)

    # Replace the Ensembl IDs with gene symbols
    gene_symbols <- gene_symbols[ensembl_ids]

    return(gene_symbols)
}

file_paths <- c(
    "/Users/jackpeplinski/CycleMix/benchmarkData/7a5c742b-d12c-4f4c-ad1d-e55649f75f7c.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/84f3485a-e4b3-49c0-8279-65762e01e0f6.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/231d025d-6b31-40da-aa38-cf618d53b544.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/e8ad2b36-b736-4ee9-889c-03555cd50165.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/fca7727d-59b3-4a5f-afa7-4d73ea824444.rds"
)

for (file_path in file_paths) {

}


simpsonIndexSeurat <- function() {
    cell_type <- seurat_cy@meta.data$cell_type
    phase <- seurat_cy$Phase
    cell_type_and_phase_df <- data.frame(cell_type, phase)
    cell_type_and_phase_table <- table(cell_type_and_phase_df)
    cell_type_and_phase_percent <- prop.table(cell_type_and_phase_table, 1)
    cell_type_and_phase_squared_percentages <- cell_type_and_phase_percent^2
    sum_of_squares_by_cell_type <- rowSums(cell_type_and_phase_squared_percentages)
    simpson_indices <- 1 - sum_of_squares_by_cell_type
}

simpsonIndexSCE <- function(sce_data, output) {
    cell_type <- colData(sce_data)$cell_type
    phase <- output$phase
    cell_type_and_phase_df <- data.frame(cell_type, phase)
    cell_type_and_phase_table <- table(cell_type_and_phase_df)
    cell_type_and_phase_percent <- prop.table(cell_type_and_phase_table, 1)
    cell_type_and_phase_squared_percentages <- cell_type_and_phase_percent^2
    sum_of_squares_by_cell_type <- rowSums(cell_type_and_phase_squared_percentages)
    simpson_indices <- 1 - sum_of_squares_by_cell_type
    return(simpson_indices)
}
