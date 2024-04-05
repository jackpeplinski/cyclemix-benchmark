library(SingleCellExperiment)
library(dplyr)
options(max.print = 20)

get_data <- function() {
    files <- list.files(
        path =
            "./benchmarkData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt",
        full.names = TRUE
    )

    # start at file 4 because the first 3 files do not have phases from ColDataXML.rds
    data <- read.table(files[4], header = TRUE, row.names = 1)
    data <- data["fpkm"]
    file_name_without_extension <- sub(".*/(GSM[0-9]+).*", "\\1", files[4])
    data <- rename(data, !!file_name_without_extension := "fpkm")

    col_data <- readRDS("./benchmarkData/ColDataXML.rds")
    for (file in files[-4]) {
        if (any(sapply(rownames(col_data), function(x) grepl(x, file)))) {
            df <- read.table(file, header = TRUE, row.names = 1)
            df <- df["fpkm"]
            file_name_without_extension <- sub(".*/(GSM[0-9]+).*", "\\1", file)
            df <- rename(df, !!file_name_without_extension := "fpkm")
            data <- cbind(data, df)
        }
    }
    return(as.matrix(data))
}

get_feature_symbol <- function() {
    # read any file to get the feature_symbol data
    data <- read.table("./benchmarkData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt", header = TRUE, row.names = 1)
    feature_symbol <- data["gene.symbol"]
    feature_symbol_df <- DataFrame(feature_symbol = feature_symbol$gene.symbol, row.names = rownames(data))
    return(feature_symbol_df)
}

get_counts <- function(data) {
    f <- function(x) if (is.numeric(x)) round(x * 100) else x
    counts <- apply(counts, c(1, 2), f)
    return(counts)
}

get_logcounts <- function(data) {
    libsizes <- colSums(data)
    size.factors <- libsizes / mean(libsizes)
    logcounts <- log2(t(t(data)) + 1)
    return(logcounts)
}

get_sce <- function() {
    data <- get_data()
    logcounts <- get_logcounts(data)
    counts <- get_counts(data)
    col_data <- readRDS("./benchmarkData/ColDataXML.rds")
    feature_symbol <- get_feature_symbol()
    sce <- SingleCellExperiment(
        assays = list(counts = counts, logcounts = logcounts),
        colData = col_data,
        rowData = feature_symbol
    )
    colnames(sce) <- colnames(counts)
    rownames(sce) <- rownames(counts)
    return(sce)
}

sce <- get_sce()


# Create a SingleCellExperiment object

# counts <- assay(sce, "counts")
# libsizes <- colSums(logcounts)
# size.factors <- libsizes / mean(libsizes)
# logcounts(sce) <- log2(t(t(logcounts)) + 1)
