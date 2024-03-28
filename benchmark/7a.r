require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("mclust")
source("R/Analysis.R") # load after cyclemix
require("stringr")
library("Seurat")
library("biomaRt")
library("tidyverse")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
library(ggplot2)
options(max.print = 20)

convert_to_ensembl <- function(gene_symbols) {
    # Convert gene symbols to Ensembl IDs
    ensembl_ids_df <- suppressMessages(suppressWarnings(
        AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = "ENSEMBL", keytype = "SYMBOL")
    ))

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
    gene_symbols_df <- suppressMessages(suppressWarnings(
        AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids, columns = "SYMBOL", keytype = "ENSEMBL")
    ))

    # If there are multiple SYMBOLs for a single ENSEMBL, keep the first one
    gene_symbols_df <- gene_symbols_df[!duplicated(gene_symbols_df$ENSEMBL), ]

    # Create a named vector of SYMBOLs with ENSEMBLs as names
    gene_symbols <- setNames(gene_symbols_df$SYMBOL, gene_symbols_df$ENSEMBL)

    # Replace the Ensembl IDs with gene symbols
    gene_symbols <- gene_symbols[ensembl_ids]

    return(gene_symbols)
}

get_cell_type_and_phase_percent <- function(cell_type, phase) {
    cell_type_and_phase_df <- data.frame(cell_type, phase)
    cell_type_and_phase_table <- table(cell_type_and_phase_df)
    cell_type_and_phase_percent <- prop.table(cell_type_and_phase_table, 1)
    return(cell_type_and_phase_percent)
}

get_simpson_index <- function(cell_type_and_phase_percent) {
    cell_type_and_phase_squared_percentages <- cell_type_and_phase_percent^2
    sum_of_squares_by_cell_type <- rowSums(cell_type_and_phase_squared_percentages)
    simpson_indices <- sum_of_squares_by_cell_type
    return(simpson_indices)
}

create_graph <- function(sce_file_path) {
    seurat_data <- readRDS(sce_file_path)

    s.genes <- HGeneSets$Whitfield$Gene[HGeneSets$Whitfield$Stage == "S"]
    s.genes <- convert_to_ensembl(as.character(s.genes))
    g2m.genes <- HGeneSets$Whitfield$Gene[HGeneSets$Whitfield$Stage == "G2M"]
    g2m.genes <- convert_to_ensembl(as.character(g2m.genes))
    seurat_cy <- CellCycleScoring(seurat_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

    sce_data <- as.SingleCellExperiment(seurat_data)
    rownames(sce_data) <- convert_to_symbols(rownames(sce_data))
    sce_data <- sce_data[!is.na(rownames(sce_data)), ]
    sce_data <- sce_data[!duplicated(rownames(sce_data)), ]
    rowData_sce_data <- DataFrame(feature_symbol = factor(rownames(sce_data)))
    rowData(sce_data) <- rowData_sce_data
    subsettedHGeneSets <- HGeneSets$Whitfield[HGeneSets$Whitfield$Stage %in% c("S", "G2M"), ]
    output <- classifyCells(sce_data, subsettedHGeneSets)

    # Assume simpsonIndexSeurat and simpsonIndexSCE are the results from the respective functions
    # Convert them to data frames
    sce_cell_type_and_phase_percent <- get_cell_type_and_phase_percent(colData(sce_data)$cell_type, output$phase)
    seurat_cell_type_and_phase_percent <- get_cell_type_and_phase_percent(seurat_cy@meta.data$cell_type, seurat_cy$Phase)
    df <- as.data.frame(cbind(simpsonIndexSCE = get_simpson_index(sce_cell_type_and_phase_percent), simpsonIndexSeurat = get_simpson_index(seurat_cell_type_and_phase_percent)))

    # Reshape the data to long format
    df_long <- df %>%
        rownames_to_column(var = "cell_type") %>%
        pivot_longer(
            cols = c(simpsonIndexSCE, simpsonIndexSeurat),
            names_to = "source",
            values_to = "simpson"
        )

    # Create the side-by-side bar plot
    p <- ggplot(df_long, aes(x = cell_type, y = simpson, fill = source)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = "Cell Type", y = "Simpson Index", fill = "Source")
    print(p)
    ggsave(paste0("my_plot_", i, ".png"), plot = p)
}


file_paths <- c(
    "/Users/jackpeplinski/CycleMix/benchmarkData/7a5c742b-d12c-4f4c-ad1d-e55649f75f7c.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/84f3485a-e4b3-49c0-8279-65762e01e0f6.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/231d025d-6b31-40da-aa38-cf618d53b544.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/e8ad2b36-b736-4ee9-889c-03555cd50165.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/fca7727d-59b3-4a5f-afa7-4d73ea824444.rds"
)

sce_file_path <- "/Users/jackpeplinski/CycleMix/benchmarkData/7a5c742b-d12c-4f4c-ad1d-e55649f75f7c.rds"
create_graph(sce_file_path)

# for (i in seq_along(file_paths)) {
#     create_graph(file_paths[i])
# }
