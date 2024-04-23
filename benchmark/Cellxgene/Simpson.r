require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("mclust")
source("R/Analysis.R") # load after cyclemix
source("./benchmark/ConvertGenes.R")
require("stringr")
library("Seurat")
library("biomaRt")
library("tidyverse")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
library(ggplot2)
options(max.print = 20)

get_cell_type_and_phase_percent <- function(cell_type, phase) {
    cell_type_and_phase_df <- data.frame(cell_type, phase)
    cell_type_and_phase_table <- table(cell_type_and_phase_df)
    cell_type_and_phase_percent <- prop.table(cell_type_and_phase_table, 1)
    return(as.data.frame.matrix(cell_type_and_phase_percent))
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

get_seurat_output <- function(seurat_data) {
    s.genes <- HGeneSets$Whitfield$Gene[HGeneSets$Whitfield$Stage == "S"]
    s.genes <- convert_hgene_symbols_to_ensembl_ids(as.character(s.genes))
    g2m.genes <- HGeneSets$Whitfield$Gene[HGeneSets$Whitfield$Stage == "G2M"]
    g2m.genes <- convert_hgene_symbols_to_ensembl_ids(as.character(g2m.genes))
    seurat_output <- CellCycleScoring(seurat_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    return(seurat_output)
}

get_simpson_graph <- function(cyclemix_cell_type_and_phase_percent, seurat_cell_type_and_phase_percent, dataset_name) {
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
        ggtitle(paste(
            "Simpson Index of ", dataset_name, " with Seurat and CycleMix"
        )) +
        geom_bar(stat = "identity", position = "dodge") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
        ) +
        labs(x = "Cell Type", y = "Simpson Index", fill = "Source")
    return(simpson_indices_graph)
}

get_cell_type_graph <- function(cyclemix_cell_type_and_phase_percent, seurat_cell_type_and_phase_percent, dataset_name) {
    # Convert the tables to data frames
    cyclemix_cell_type_and_phase_percent$cell_type <- rownames(cyclemix_cell_type_and_phase_percent)
    cyclemix_cell_type_and_phase_percent$source <- "CycleMix"

    seurat_cell_type_and_phase_percent$cell_type <- rownames(seurat_cell_type_and_phase_percent)
    seurat_cell_type_and_phase_percent$source <- "Seurat"

    # Add missing columns to both data frames
    cyclemix_cell_type_and_phase_percent$G1 <- 0
    seurat_cell_type_and_phase_percent$None <- 0

    df_long <- rbind(cyclemix_cell_type_and_phase_percent, seurat_cell_type_and_phase_percent)

    df_long <- tidyr::pivot_longer(df_long, cols = c(G1, G2M, S, None), names_to = "phase", values_to = "percent")

    p <- ggplot(df_long) +
        ggtitle(paste(
            "Predicted Cell Types of ", dataset_name, " with Seurat and CycleMix"
        )) +
        geom_bar(aes(x = source, y = percent, fill = phase),
            position = "stack",
            stat = "identity"
        ) +
        scale_fill_manual(values = c("None" = "grey", "G1" = "darkblue", "S" = "#dddd26", "G2M" = "#ff4848")) +
        facet_grid(~cell_type, switch = "x") +
        theme(
            strip.placement = "outside",
            strip.background = element_rect(fill = NA, color = "white"),
            panel.spacing = unit(.01, "cm"),
            strip.text = element_text(size = 16, angle = 90, hjust = 1),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.title = element_text(size = 16),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
        ) +
        labs(x = "Source", y = "Percent", fill = "Phase")
    return(p)
}

save_graph <- function(file_prefix, file_path, file_name, p, width = 10, height = 10, units = "in") {
    print(paste0("Saving graph ", file_prefix, file_name, ".png"))
    ggsave(paste0(file_prefix, file_name, ".png"), plot = p, width = width, height = height, units = units)
}

process_seurat_object <- function(file_path) {
    seurat_data <- readRDS(file_path)

    seurat_output <- get_seurat_output(seurat_data)
    cyclemix_output <- get_cyclemix_output(as.SingleCellExperiment(seurat_data))

    cell_types <- seurat_output@meta.data$cell_type
    cyclemix_cell_type_and_phase_percent <- get_cell_type_and_phase_percent(cell_types, cyclemix_output$phase)
    seurat_cell_type_and_phase_percent <- get_cell_type_and_phase_percent(cell_types, seurat_output$Phase)

    file_name <- tools::file_path_sans_ext(basename(file_path))
    dataset_name <- dataset_name_map[file_name]

    simpson_indices_graph <- get_simpson_graph(cyclemix_cell_type_and_phase_percent, seurat_cell_type_and_phase_percent, dataset_name)
    save_graph("output/simpson_", file_path, file_name, simpson_indices_graph)

    cell_type_graph <- get_cell_type_graph(cyclemix_cell_type_and_phase_percent, seurat_cell_type_and_phase_percent, dataset_name)
    save_graph("output/cell_type_", file_path, file_name, cell_type_graph, width = 20, height = 10, units = "in")
}

process_seurat_objects <- function(file_paths) {
    for (i in seq_along(file_paths)) {
        file_path <- file_paths[i]
        process_seurat_object(file_path)
    }
}

# single_file("/Users/jackpeplinski/CycleMix/benchmarkData/7a5c742b-d12c-4f4c-ad1d-e55649f75f7c.rds")

dataset_name_map <- list(
    "e8ad2b36-b736-4ee9-889c-03555cd50165" = "Bondoc et al. (GSE180665)",
    "84f3485a-e4b3-49c0-8279-65762e01e0f6" = "Wu et al. (GSE176078)",
    "231d025d-6b31-40da-aa38-cf618d53b544" = "Zhang et al.",
    "7a5c742b-d12c-4f4c-ad1d-e55649f75f7c" = "Darmanis et al.",
    "fca7727d-59b3-4a5f-afa7-4d73ea824444" = "Chan et al."
)

process_seurat_objects(c(
    "./benchmarkData/7a5c742b-d12c-4f4c-ad1d-e55649f75f7c.rds",
    "./benchmarkData/84f3485a-e4b3-49c0-8279-65762e01e0f6.rds",
    "./benchmarkData/231d025d-6b31-40da-aa38-cf618d53b544.rds",
    "./benchmarkData/e8ad2b36-b736-4ee9-889c-03555cd50165.rds",
    "./benchmarkData/fca7727d-59b3-4a5f-afa7-4d73ea824444.rds"
))
