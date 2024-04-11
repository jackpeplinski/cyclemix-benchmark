require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")
library("Seurat")
library("biomaRt")
library("org.Mm.eg.db")
options(max.print = 20)

MSeuratGeneSet <- readRDS("./benchmarkData/MSeuratGeneSet.RDS")
seurat_mouse_orth <- readRDS("./benchmarkData/SeuratCC_toMmus_ortho.rds")

synthesize_gse_42268 <- function(sce) {
    # create synthetic data
    synth_g1 <- synthesize_gse_42268_by_type(sce, "G1")
    synth_g2 <- synthesize_gse_42268_by_type(sce, "G2/M")
    synth_s <- synthesize_gse_42268_by_type(sce, "S")
    g1_names <- colnames(synth_g1)
    g2_names <- colnames(synth_g2)
    s_names <- colnames(synth_s)
    combined <- cbind(synth_g1, synth_g2, synth_s)

    sceMix <- SingleCellExperiment(
        assays = list(counts = combined),
        colData = DataFrame(
            row.names = c(g1_names, g2_names, s_names),
            cell_type1 = factor(c(rep("G1", length(g1_names)), rep("G2/M", length(g2_names)), rep("S", length(
                s_names
            ))))
        ),
        rowData = rowData(sce)
    )
    libsizes <- colSums(combined)
    size.factors <- libsizes / mean(libsizes)
    logcounts(sceMix) <- log2(t(t(combined)) + 1)
    return(sceMix)
}

synthesize_gse_42268_by_type <- function(sce, cell_type) {
    cells <- colnames(sce)[colData(sce)$cell_type1 == cell_type]
    phase_df <- assays(sce)$counts[, cells]
    synth <- average_phase_df(phase_df)
    return(synth)
}

average_phase_df <- function(phase_df) {
    n <- ncol(phase_df)
    synthetic_cells <- data.frame(matrix(NA, nrow = nrow(phase_df), ncol = 0))

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            col1 <- phase_df[, i]
            col2 <- phase_df[, j]
            col_name <- paste(colnames(phase_df)[i], "-", colnames(phase_df)[j], sep = "")
            synthetic_cells[[col_name]] <- rowMeans(cbind(col1, col2))
        }
    }

    return(cbind(synthetic_cells, phase_df))
}

# Function to calculate the average of two adjacent columns
average_of_adjacent <- function(df) {
    # Determine the number of columns
    num_cols <- ncol(df)

    # Check if the number of columns is odd, and if so, remove the last column
    if (num_cols %% 2 != 0) {
        num_cols <- num_cols - 1
        df <- df[, -num_cols]
    }

    new_df <- df
    for (column_i in 1:num_cols) {
        # add columns by putting square brackets with the column name
        average <- as.data.frame(rowMeans(s_df[, c(column_i, column_i + 1)]))
        new_df <- rbind(new_df, average)
    }

    # Calculate the average of two adjacent columns and store in a new dataframe
    new_df <- as.data.frame(matrix(rowMeans(matrix(df, ncol = 2, byrow = TRUE)), ncol = (num_cols / 2)))

    # Generate a header for the new column (e.g., "Avg_GSM1_GSM2")
    col_header <- paste("Avg", colnames(df)[seq(1, num_cols, 2)], colnames(df)[seq(2, num_cols, 2)], sep = "_")

    # Set the new column header
    colnames(new_df) <- col_header

    return(new_df)
}

validate_gse_42268 <- function() {
    sce <- get_sce()
    invalid_row_counts <- c()
    for (i in colnames(logcounts(sce))) { # names of files
        # print(colnames(logcounts(sce)))
        invalid_row_count <- 0
        # if (i == "GSM1036501") {
        #     print(logcounts(sce)[i])
        # }
        file_log_counts <- logcounts(sce)[i]
        # print(colnames(fileLogCounts))
        colnames(file_log_counts)

        logcounts_rows <- row.names(file_log_counts)
        cyclone_rows <- row.names(MGeneSets$Cyclone)
        for (cycloneRow in cyclone_rows) {
            if (is.element(cycloneRow, logcounts_rows)) {
                # print(head(fileLogCounts))
                if (file_log_counts[cycloneRow, ] < 5) {
                    invalid_row_count <- invalid_row_count + 1
                    # print(invalidRowCount)
                }
                # print(cycloneRow)
                # print(type(fileLogCounts[fileLogCounts[[1]] > 5, ]))
                # print(nrow(fileLogCounts))
            }
        }
        invalid_row_counts <- append(invalid_row_counts, invalid_row_count
        / nrow(MGeneSets$Cyclone) * 100)
        # print(invalidRowCounts)
    }
    return(invalid_row_counts)
    # try subsetting MGeneSets$Cyclone by removing rows with Dir of -1 => this did improve it... but not a lot
    # look at average logcounts, the rows that have genes from MGeneSets$Cyclone. If the values are small (less than 5) they are not good indicators of the phase and how many are there. => 14 cells had poor logcounts, but 46 were not classified.
}
