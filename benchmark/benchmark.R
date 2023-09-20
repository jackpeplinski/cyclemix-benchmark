require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")
library("Seurat")
library("biomaRt")
library("org.Mm.eg.db")

"
GSE42268 Dataset Notes:
- Don't need to log normalize values
- Each file is a cell
- Round fpkm for counts
- Put fpkm into the logcounts
"
format_data <- function() {
    data <- readRDS("./benchmarkData/SeuratCC_toMmus_ortho.rds")

    phase <- c("S")

    # Determine the number of rows in the data frame
    num_rows <- length(data$mmus_s)

    # Repeat the phase values to match the number of rows
    phase <- rep_len(phase, num_rows)

    # Combine the columns using cbind()
    result <- cbind(data$mmus_s, phase)

    # Convert the result to a data frame
    result <- as.data.frame(result)

    # Assign column names
    colnames(result) <- c("Gene", "Stage")

    ##### Phase 2
    phase <- c("G2M")

    # Determine the number of rows in the data frame
    num_rows <- length(data$mmus_g2m)

    # Repeat the phase values to match the number of rows
    phase <- rep_len(phase, num_rows)

    # Combine the columns using cbind()
    result2 <- cbind(data$mmus_g2m, phase)

    # Convert the result to a data frame
    result2 <- as.data.frame(result2)

    # Assign column names
    colnames(result2) <- c("Gene", "Stage")

    results <- rbind(result, result2)
    results$Gene <- results$Gene
    results$Dir <- 1
    results$Gene <- as.factor(results$Gene)
    results$Stage <- as.factor(results$Stage)

    ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

    ensembl_ids <- getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "external_gene_name",
        values = results$Gene,
        mart = ensembl
    )
    ensembl_ids$external_gene_name <- ensembl_ids$external_gene_name

    merged_df <- merge(results, ensembl_ids, by.x = "Gene", by.y = "external_gene_name", all.x = TRUE)

    rownames(merged_df) <- merged_df$ensembl_gene_id

    merged_df <- subset(merged_df, select = -ensembl_gene_id)
    return(merged_df)
}
# MSeuratGeneSet <- format_data()
# seurat_mouse_orth <- readRDS("./benchmarkData/SeuratCC_toMmus_ortho.rds")

format_gse_42268 <- function() {
    # counts should have a structure like:
    # row.name                 gsmcellname
    #                          [fpkm]

    # load counts info; doesn't matter which file, could be any
    counts <- read.table(
        str_interp("./benchmarkData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )

    # adds row.name with ensembleID
    rownames(counts) <- counts$id

    # adds gene.symbol with gene symbol
    row_data <- DataFrame(
        row.names = counts$id,
        feature_symbol = factor(counts$gene.symbol)
    )

    # remove unneeded info. counts is empty with only rownames = ensembl id now
    counts <- subset(counts,
        select = -c(
            id,
            gene.symbol,
            fpkm
        )
    )

    # duplicate counts variable for logcounts
    logcounts <- counts

    # copy/pasted this from gseParse
    col_data_xml <- data.frame(
        cell_type1 = c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1"),
        gsm = c("GSM1036483", "GSM1036484", "GSM1036485", "GSM1036486", "GSM1036487", "GSM1036488", "GSM1036489", "GSM1036490", "GSM1036491", "GSM1036492", "GSM1036493", "GSM1036494", "GSM1036495", "GSM1036496", "GSM1036500", "GSM1036501", "GSM1036502", "GSM1036503", "GSM1036504", "GSM1036505", "GSM1036506", "GSM1036507", "GSM1036508", "GSM1036509", "GSM1036510", "GSM1036511", "GSM1036512", "GSM1036513", "GSM1036522", "GSM1036523", "GSM1036524", "GSM1036525", "GSM1036526", "GSM1036527", "GSM1036528", "GSM1036529", "GSM1036530", "GSM1036531", "GSM1036532", "GSM1036533", "GSM1036534", "GSM1036535", "GSM1036536", "GSM1036537", "GSM1036538", "GSM1036539", "GSM1036540", "GSM1036541", "GSM1036542", "GSM1036543", "GSM1036544", "GSM1036545", "GSM1036546", "GSM1036547", "GSM1036548", "GSM1036549", "GSM1036550", "GSM1036551", "GSM1036552", "GSM1036553", "GSM1036554", "GSM1036555", "GSM1036556")
    )

    files <- list.files(
        path =
            "./benchmarkData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    # for each file
    for (file in files) {
        # if the file name matches the col_data
        gsm <- substr(file, 1, 10)
        if (is.element(gsm, col_data_xml$gsm)) {
            # get the data for the file
            counts_data <- read.table(
                str_interp("./benchmarkData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )

            # format the data
            format_counts <- function(count_type) {
                # round data if counts
                if (count_type == "counts") {
                    data <- data.frame(
                        lapply(
                            counts_data,
                            function(x) if (is.numeric(x)) round(x) else x
                        )
                    )
                } else {
                    data <- counts_data
                }

                # sets rownames to ids
                rownames(data) <- data$id

                # remove uneeded data
                data <- subset(data,
                    select = -c(
                        id,
                        gene.symbol
                    )
                )

                colnames(data)[1] <- gsm
                return(data)
            }

            # the values for the gene expression
            logcounts_data <- format_counts("logcounts")
            counts_data <- format_counts("counts")

            # the phases of cells
            col_data <- DataFrame(
                row.names = col_data_xml$gsm,
                cell_type1 = factor(col_data_xml$cell_type1)
            )
            # verify that the correct phase is matched
            # print(col_data["GSM1036531", ])

            # combine phases and expressions
            logcounts <- cbind(logcounts, logcounts_data)
            counts <- cbind(counts, counts_data)
        }
    }

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data,
        rowData = row_data
    )
    counts <- assay(sce, "counts")
    libsizes <- colSums(logcounts)
    size.factors <- libsizes / mean(libsizes)
    logcounts(sce) <- log2(t(t(logcounts)) + 1)

    # because this for some reason isn't working with ensembleids, remove the duplicates and use gene idsa
    df <- rowData(sce)
    duplicated_rows <- which(duplicated(df$feature_symbol))
    row_names <- df[-duplicated_rows, ]
    counts <- assays(sce)$counts[-duplicated_rows, ]
    row.names(counts) <- row_names
    logcounts <- assays(sce)$logcounts[-duplicated_rows, ]
    row.names(logcounts) <- row_names
    row_data <- DataFrame(
        row.names = row_names,
        feature_symbol = factor(row_names)
    )
    col_data <- as.data.frame(colData(sce))[, , FALSE]
    col_data <- DataFrame(
        row.names = rownames(col_data),
        cell_type1 = factor(col_data$cell_type1)
    )

    sce <- SingleCellExperiment(
        assays = list(counts = counts, logcounts = logcounts),
        colData = col_data,
        rowData = row_data
    )

    # create synthetic data
    synth_g1 <- synthesize_gse_42268(sce, "G1")
    synth_g2 <- synthesize_gse_42268(sce, "G2/M")
    synth_s <- synthesize_gse_42268(sce, "S")
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
        rowData = row_data
    )
    libsizes <- colSums(combined)
    size.factors <- libsizes / mean(libsizes)
    logcounts(sceMix) <- log2(t(t(combined)) + 1)
    return(sceMix)
}

synthesize_gse_42268 <- function(sce, cell_type) {
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

classify_gse_42268 <- function() {
    gse_sce <<- format_gse_42268()
    cat("===GSE 42268 | CycleMix | MGeneSets$Cyclone===\n")
    gse_cm_cy <<- classifyCells(gse_sce, MGeneSets$Cyclone)
    print(table(factor(gse_cm_cy$phase), gse_sce$cell_type1))
    cat("======GSE 42268 | CycleMix | MSeuratGeneSet===\n")
    gse_cm_se <<- classifyCells(gse_sce, MSeuratGeneSet)
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
# classify_gse_42268()

validate_gse_42268 <- function() {
    sce <- format_gse_42268()
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

    # get cell phases from python script
    col_data_xml <- data.frame(
        cell_type1 = c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1"),
        gsm = c("GSM1036483", "GSM1036484", "GSM1036485", "GSM1036486", "GSM1036487", "GSM1036488", "GSM1036489", "GSM1036490", "GSM1036491", "GSM1036492", "GSM1036493", "GSM1036494", "GSM1036495", "GSM1036496", "GSM1036500", "GSM1036501", "GSM1036502", "GSM1036503", "GSM1036504", "GSM1036505", "GSM1036506", "GSM1036507", "GSM1036508", "GSM1036509", "GSM1036510", "GSM1036511", "GSM1036512", "GSM1036513", "GSM1036522", "GSM1036523", "GSM1036524", "GSM1036525", "GSM1036526", "GSM1036527", "GSM1036528", "GSM1036529", "GSM1036530", "GSM1036531", "GSM1036532", "GSM1036533", "GSM1036534", "GSM1036535", "GSM1036536", "GSM1036537", "GSM1036538", "GSM1036539", "GSM1036540", "GSM1036541", "GSM1036542", "GSM1036543", "GSM1036544", "GSM1036545", "GSM1036546", "GSM1036547", "GSM1036548", "GSM1036549", "GSM1036550", "GSM1036551", "GSM1036552", "GSM1036553", "GSM1036554", "GSM1036555", "GSM1036556")
    )

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
    col_data_xml <- data.frame(
        cell_type1 = c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1"),
        gsm = c("GSM1036483", "GSM1036484", "GSM1036485", "GSM1036486", "GSM1036487", "GSM1036488", "GSM1036489", "GSM1036490", "GSM1036491", "GSM1036492", "GSM1036493", "GSM1036494", "GSM1036495", "GSM1036496", "GSM1036500", "GSM1036501", "GSM1036502", "GSM1036503", "GSM1036504", "GSM1036505", "GSM1036506", "GSM1036507", "GSM1036508", "GSM1036509", "GSM1036510", "GSM1036511", "GSM1036512", "GSM1036513", "GSM1036522", "GSM1036523", "GSM1036524", "GSM1036525", "GSM1036526", "GSM1036527", "GSM1036528", "GSM1036529", "GSM1036530", "GSM1036531", "GSM1036532", "GSM1036533", "GSM1036534", "GSM1036535", "GSM1036536", "GSM1036537", "GSM1036538", "GSM1036539", "GSM1036540", "GSM1036541", "GSM1036542", "GSM1036543", "GSM1036544", "GSM1036545", "GSM1036546", "GSM1036547", "GSM1036548", "GSM1036549", "GSM1036550", "GSM1036551", "GSM1036552", "GSM1036553", "GSM1036554", "GSM1036555", "GSM1036556")
    )

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

format_emtab_2805 <- function() {
    emtab_2805_file <- function(file_name) {
        # get values
        counts <- read.table(
            str_interp("./benchmarkData/E-MTAB-2805/E-MTAB-2805.processed.1/${file_name}.txt"),
            header = TRUE
        )
        counts$AssociatedGeneName <- counts$AssociatedGeneName
        counts <- counts[!duplicated(counts$AssociatedGeneName), ]
        counts <- na.omit(counts)

        # set rownames
        rownames(counts) <- counts$AssociatedGeneName

        # remove unused values
        counts <- head(counts, -97)

        # build needed df for sce
        row_data <- DataFrame(
            row.names = counts$AssociatedGeneName,
            feature_symbol = factor(counts$AssociatedGeneName)
        )

        # build needed df for sce
        col_data <- DataFrame(
            row.names = colnames(counts)[5:100],
            Species = factor("Mus musculus"),
            cell_type1 = factor(
                substr(colnames(counts)[5:100], 1, 2)
            ),
            Source = factor("ESC")
        )

        # remove unneeded data
        counts <- subset(counts,
            select = -c(
                EnsemblGeneID,
                EnsemblTranscriptID,
                AssociatedGeneName,
                GeneLength
            )
        )

        # build sce
        sce <- SingleCellExperiment(
            assays = list(counts = counts),
            colData = col_data,
            rowData = row_data
        )
        counts <- assay(sce, "counts")
        libsizes <- colSums(counts)
        size.factors <- libsizes / mean(libsizes)
        logcounts(sce) <- log2(t(t(counts) / size.factors) + 1)
        return(sce)
    }

    emtab_2805_res <- cbind(
        emtab_2805_file("G1_singlecells_counts"),
        emtab_2805_file("G2M_singlecells_counts"),
        emtab_2805_file("S_singlecells_counts")
    )
    return(emtab_2805_res)
}

classify_emtab_2805 <- function() {
    cat("===EMTAB 2805 | CycleMix | MGeneSets$Cyclone===\n")
    emtab_sce <<- format_emtab_2805()
    emtab_cm_cy <<- classifyCells(emtab_sce, MGeneSets$Cyclone)
    print(table(factor(emtab_cm_cy$phase), emtab_sce$cell_type1))
    cat("===EMTAB 2805 | CycleMix | MSeuratGeneSet===\n")
    emtab_cm_se <<- classifyCells(emtab_sce, MSeuratGeneSet)
    emtab_cm_se_table <- table(factor(emtab_cm_se$phase), emtab_sce$cell_type1)
    print(emtab_cm_se_table)

    emtab_seurat <- as.Seurat(emtab_sce)
    emtab_seurat <- NormalizeData(emtab_seurat)
    emtab_seurat <- FindVariableFeatures(emtab_seurat, selection.method = "vst")
    emtab_seurat <- ScaleData(emtab_seurat, features = rownames(emtab_seurat))
    emtab_seurat <<- RunPCA(emtab_seurat, features = VariableFeatures(emtab_seurat), ndims.print = 6:10, nfeatures.print = 10)
    cat("===EMTAB 2805 | Seurat | MGeneSets$Cyclone===\n")
    s.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "S"]
    g2m.genes <- MGeneSets$Cyclone$Gene[MGeneSets$Cyclone$Stage == "G2M"]
    emtab_seurat_cy <<- CellCycleScoring(emtab_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(emtab_seurat_cy[[]]$Phase, emtab_seurat[[]]$orig.ident))
    cat("===EMTAB 2805 | Seurat | MSeuratGeneSet===\n")
    s.genes <- seurat_mouse_orth$mmus_s
    g2m.genes <- seurat_mouse_orth$mmus_g2m
    emtab_seurat_se <<- CellCycleScoring(emtab_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    print(table(emtab_seurat_se[[]]$Phase, emtab_seurat[[]]$orig.ident))
}
# classify_emtab_2805()
