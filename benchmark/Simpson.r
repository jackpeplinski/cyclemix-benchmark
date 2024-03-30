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

convert_gene_symbols_to_ensembl_ids <- function(gene_symbols) {
    ensembl_ids_df <- suppressMessages(suppressWarnings(
        AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = "ENSEMBL", keytype = "SYMBOL")
    ))

    ensembl_ids_df <- ensembl_ids_df[!is.na(ensembl_ids_df$ENSEMBL), ]

    # If there are multiple ENSEMBLs for a single SYMBOL, keep the first one
    ensembl_ids_df <- ensembl_ids_df[!duplicated(ensembl_ids_df$SYMBOL), ]

    return(ensembl_ids_df$ENSEMBL)
}

convert_ensemble_ids_to_gene_symbols <- function(ensembl_ids) {
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

get_cyclemix_output <- function(sce_data) {
    rownames(sce_data) <- convert_ensemble_ids_to_gene_symbols(rownames(sce_data))
    sce_data <- sce_data[!is.na(rownames(sce_data)), ]
    sce_data <- sce_data[!duplicated(rownames(sce_data)), ]
    rowData_sce_data <- DataFrame(feature_symbol = factor(rownames(sce_data)))
    rowData(sce_data) <- rowData_sce_data
    subsettedHGeneSets <- HGeneSets$Whitfield[HGeneSets$Whitfield$Stage %in% c("S", "G2M"), ]
    cyclemix_output <- classifyCells(sce_data, subsettedHGeneSets)
    return(cyclemix_output)
}



create_simpson_and_cell_type_graph <- function(sce_file_path) {
    seurat_data <- readRDS(sce_file_path)

    # Run seurat on the data
    s.genes <- HGeneSets$Whitfield$Gene[HGeneSets$Whitfield$Stage == "S"]
    s.genes <- convert_gene_symbols_to_ensembl_ids(as.character(s.genes))
    g2m.genes <- HGeneSets$Whitfield$Gene[HGeneSets$Whitfield$Stage == "G2M"]
    g2m.genes <- convert_gene_symbols_to_ensembl_ids(as.character(g2m.genes))
    seurat_output <- CellCycleScoring(seurat_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    sce_data <- as.SingleCellExperiment(seurat_data)
    cyclemix_output <- get_cyclemix_output(sce_data)

    cyclemix_cell_type_and_phase_percent <- get_cell_type_and_phase_percent(colData(sce_data)$cell_type, cyclemix_output$phase)
    seurat_cell_type_and_phase_percent <- get_cell_type_and_phase_percent(seurat_output@meta.data$cell_type, seurat_output$Phase)
    simpson_indices_by_cell_type <- as.data.frame(cbind(CycleMix = get_simpson_index(cyclemix_cell_type_and_phase_percent), Seurat = get_simpson_index(seurat_cell_type_and_phase_percent)))

    # Reshape the data to long format
    simpson_indices_by_cell_type_long <- simpson_indices_by_cell_type %>%
        rownames_to_column(var = "cell_type") %>%
        pivot_longer(
            cols = c(CycleMix, Seurat),
            names_to = "source",
            values_to = "simpson"
        )

    simpson_indices_graph <- ggplot(simpson_indices_by_cell_type_long, aes(x = cell_type, y = simpson, fill = source)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = "Cell Type", y = "Simpson Index", fill = "Source")
    datafile_name <- tools::file_path_sans_ext(basename(sce_file_path))
    ggsave(paste0("output/simpson_", datafile_name, ".png"), plot = simpson_indices_graph)

    # Convert the tables to data frames
    sce_df <- as.data.frame.matrix(cyclemix_cell_type_and_phase_percent)
    sce_df$cell_type <- rownames(sce_df)
    sce_df$source <- "sce"

    seurat_df <- as.data.frame.matrix(seurat_cell_type_and_phase_percent)
    seurat_df$cell_type <- rownames(seurat_df)
    seurat_df$source <- "seurat"

    # Add missing columns to both data frames
    sce_df$G1 <- 0
    seurat_df$None <- 0

    df_long <- rbind(sce_df, seurat_df)

    df_long <- tidyr::pivot_longer(df_long, cols = c(G1, G2M, S, None), names_to = "phase", values_to = "percent")

    p <- ggplot(df_long) +
        geom_bar(aes(x = source, y = percent, fill = phase),
            position = "stack",
            stat = "identity"
        ) +
        facet_grid(~cell_type, switch = "x") +
        theme(
            strip.placement = "outside",
            strip.background = element_rect(fill = NA, color = "white"),
            panel.spacing = unit(.01, "cm")
        )
    ggsave(paste0("output/cell_type_", datafile_name, ".png"), plot = p, width = 20, height = 10, units = "in")
}




file_paths <- c(
    "/Users/jackpeplinski/CycleMix/benchmarkData/7a5c742b-d12c-4f4c-ad1d-e55649f75f7c.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/84f3485a-e4b3-49c0-8279-65762e01e0f6.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/231d025d-6b31-40da-aa38-cf618d53b544.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/e8ad2b36-b736-4ee9-889c-03555cd50165.rds",
    "/Users/jackpeplinski/CycleMix/benchmarkData/fca7727d-59b3-4a5f-afa7-4d73ea824444.rds"
)

sce_file_path <- "/Users/jackpeplinski/CycleMix/benchmarkData/7a5c742b-d12c-4f4c-ad1d-e55649f75f7c.rds"
create_simpson_and_cell_type_graph(sce_file_path)

# for (i in seq_along(file_paths)) {
#     create_simpson_and_cell_type_graph(file_paths[i])
# }
