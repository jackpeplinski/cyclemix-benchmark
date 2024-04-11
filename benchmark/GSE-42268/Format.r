library(SingleCellExperiment)
library(dplyr)
library("org.Mm.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
options(max.print = 20)

does_file_have_phase <- function(col_data, file) {
    return(any(sapply(rownames(col_data), function(x) grepl(x, file))))
}

# for example, given GSM1036480_EB5K_01.txt return GSM1036480
get_gsm_number <- function(file) {
    return(sub(".*/(GSM[0-9]+).*", "\\1", file))
}

get_data_from_file <- function(file) {
    data <- read.table(file, header = TRUE, row.names = 1)
    data <- data["fpkm"]
    gsm_number <- get_gsm_number(file)
    data <- rename(data, !!gsm_number := "fpkm")
    return(data)
}

get_data_from_files <- function() {
    files <- list.files(
        path =
            "./benchmarkData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt",
        full.names = TRUE
    )

    # start at file 4 because the first 3 files do not have phases from ColDataXML.rds
    data <- get_data_from_file(files[4])

    col_data <- readRDS("./benchmarkData/ColDataXML.rds")
    for (file in files[-4]) {
        if (does_file_have_phase(col_data, file)) {
            df <- get_data_from_file(file)
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
    data <- apply(data, c(1, 2), f)
    return(data)
}

get_logcounts <- function(data) {
    libsizes <- colSums(data)
    size.factors <- libsizes / mean(libsizes)
    logcounts <- log2(t(t(data)) + 1)
    return(logcounts)
}

get_sce <- function() {
    data <- get_data_from_files()
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
