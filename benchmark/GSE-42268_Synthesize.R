require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")
library("Seurat")
library("biomaRt")
library("org.Mm.eg.db")
options(max.print = 20)

"
GSE42268 Dataset Notes:
- Don't need to log normalize values
- Each file is a cell
- Round fpkm for counts
- Put fpkm into the logcounts
"
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

classify_gse_42268 <- function(gse_sce) {
    cat("===GSE 42268 | CycleMix | MGeneSets$Cyclone===\n")
    gse_cm_cy <<- classifyCells(gse_sce, MGeneSets$Cyclone)
    # plotMixture(gse_cm_cy$fit[["G1"]], BIC = TRUE)
    # plotMixture(gse_cm_cy$fit[["G2M"]], BIC = TRUE)
    # plotMixture(gse_cm_cy$fit[["S"]], BIC = TRUE)
    print(table(factor(gse_cm_cy$phase), gse_sce$cell_type1))
    cat("======GSE 42268 | CycleMix | MSeuratGeneSet===\n")
    gse_cm_se <<- classifyCells(gse_sce, MSeuratGeneSet)
    # plotMixture(gse_cm_se$fit[["G1"]], BIC = TRUE)
    # plotMixture(gse_cm_se$fit[["G2M"]], BIC = TRUE)
    # plotMixture(gse_cm_se$fit[["S"]], BIC = TRUE)
    print(table(factor(gse_cm_se$phase), gse_sce$cell_type1))
    gse_seurat <- as.Seurat(gse_sce)
    gse_seurat <- NormalizeData(gse_seurat)
    gse_seurat <- FindVariableFeatures(gse_seurat, selection.method = "vst")
    gse_seurat <- ScaleData(gse_seurat, features = rownames(gse_seurat))
    gse_seurat <<- RunPCA(gse_seurat, features = VariableFeatures(gse_seurat), ndims.print = 6:10, nfeatures.print = 10)
    cat("===GSE 42268 | Seurat | MGeneSets$Cyclone===\n")
    s.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "S"]
    g2m.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "G2M"]
    gse_seurat_cy <<- CellCycleScoring(gse_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(gse_seurat_cy[[]]$Phase, gse_seurat[[]]$cell_type1))
    cat("======GSE 42268 | Seurat | MSeuratGeneSet===\n")
    s.genes <- seurat_mouse_orth$mmus_s
    g2m.genes <- seurat_mouse_orth$mmus_g2m
    gse_seurat_se <<- CellCycleScoring(gse_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(gse_seurat_se[[]]$Phase, gse_seurat[[]]$cell_type1))
}
gse_sce <<- get_sce()
cat("***Non-synthetic***\n")
classify_gse_42268(gse_sce)
cat("***Synthetic***\n")
classify_gse_42268(synthesize_gse_42268(gse_sce))

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

# for each gene in the cyclemix
# wilcox.test(#vector of all the values of the cells that are G2M, #same but for S/G1)
wilcox_gse_42268 <- function() {
    "
    get row names from file
    for each ensemble id
        for each file
            get the value of that row name
            if that file is G2M store that value in a vector
            else store that value in another vector
    run a wilcox test
    "

    # get ensemble_ids
    counts <- read.table(
        str_interp("./benchmarkData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )
    ensemble_ids <- counts$id

    # get list of files
    files <- list.files(
        path =
            "./benchmarkData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    col_data_xml <- readRDS("./benchmarkData/ColDataXML.rds")


    # for each id
    for (ensemble_id in ensemble_ids) {
        # if the element is in Cyclone set
        if (is.element(ensemble_id, rownames(MGeneSets$Cyclone))) {
            g2m <- c()
            g1 <- c()
            s <- c()

            # for each file
            for (file in files) {
                # get the file name
                gsm <- substr(file, 1, 10)

                # get values for the file
                counts <- read.table(
                    str_interp("./benchmarkData/GSE42268/GSE42268_RAW/${file}"),
                    header = TRUE
                )

                # get the values for that ID
                fpkm <- counts[counts$id == ensemble_id, ]$fpkm

                # print(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1)
                # print(length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0)

                # add the fpkm value for the ensemble for that cell
                if (length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0 &&
                    col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1 == "G2/M") {
                    g2m <- append(g2m, fpkm)
                } else if (length(col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1) > 0 &&
                    col_data_xml[col_data_xml$gsm == gsm, ]$cell_type1 == "S") {
                    g2m <- append(s, fpkm)
                } else {
                    g1 <- append(g1, fpkm)
                }
            }
            print(g2m)
            print(paste(ensemble_id, wilcox.test(g2m, s)$p.value))
            print(paste(ensemble_id, wilcox.test(g2m, g1)$p.value))
        }
    }
}

wilcox_fast_gse_42268 <- function() {
    "
        row.name (ensembleID)    gsmcellname
        ENSMUSG00000000049       0.002000
    "

    # set row names to ensemble ids
    counts <- read.table(
        str_interp("./benchmarkData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )
    ensemble_ids <- counts$id
    df <- data.frame(row.names = ensemble_ids)

    # set columns with fpkms
    files <- list.files(
        path =
            "./benchmarkData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    "
        gsm          cell_type1
        GSM1036483   G1
        ...
    "
    # not all cells have phases
    col_data_xml <- readRDS("./benchmarkData/ColDataXML.rds")
    gsms <- col_data_xml$gsm
    for (file in files) {
        gsm <- substr(file, 1, 10)
        if (is.element(gsm, gsms)) {
            counts_data <- read.table(
                str_interp("./benchmarkData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )
            counts_data <- subset(counts_data,
                select = -c(
                    id,
                    gene.symbol
                )
            )
            df[, gsm] <- counts_data
        }
    }

    print(head(df))

    # wilcox.test(#vector of all the values of the cells for each gene that are G2M, #same but for S/G1)
    # g2m <- c()
    # g1 <- c()
    # s <- c()
    gsmsG2M <- col_data_xml[col_data_xml$cell_type1 == "G2/M", ]$gsm
    gsmsG1 <- col_data_xml[col_data_xml$cell_type1 == "G1", ]$gsm
    gsmsS <- col_data_xml[col_data_xml$cell_type1 == "S", ]$gsm

    countsG2M <- df[gsmsG2M]
    countsG1 <- df[gsmsG1]
    countsS <- df[gsmsS]

    print("Writing to file...")
    file <- file("output.txt", "w")
    for (row in rownames(countsG2M)) {
        x <- countsG2M[row, ]
        y <- countsG1[row, ]
        x <- 5
        # also display the mean and output to a file so that we can sort it by pvalue
        pval <- p.adjust(wilcox.test(as.numeric(unlist((x))), as.numeric(unlist(y)))$p.value, method = "fdr")
        line <- paste(row, pval, sep = " ")
        # print(line)
        writeLines(line, file)
        # file = file("wilcox.txt")
    }
    close(file)
    print("Writing done.")
    # adjusted p-value
    # how many are below 0.05

    # print(head(countsG2M))
    # print(wilcox.test(as.numeric(unlist(countsG2M)), as.numeric(unlist(countsG1)))$p.value)
    # print(wilcox.test(as.numeric(unlist(countsG2M)), as.numeric(unlist(countsS)))$p.value)

    readFile <- function() {
        df <- read.table(
            "output.txt",
            header = TRUE
        )
        return(df)
    }
    # df <- readFile() # nrow(df[df$p_value<0.05,]), mean(df$p_value)
    # df2 <- df[order(df$p_value), ]
}
