require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")

e_mtab <- function(file_name) {
    counts <- read.table(
        str_interp("./jackData/E-MTAB-2805.processed.1/${file_name}.txt"),
        header = TRUE
    )

    rownames(counts) <- counts$EnsemblGeneID

    counts <- head(counts, -97)
    row_data <- DataFrame(
        row.names = counts$EnsemblGeneID,
        feature_symbol = factor(counts$AssociatedGeneName)
    )

    col_data <- DataFrame(
        row.names = colnames(counts)[5:100],
        Species = factor("Mus musculus"),
        cell_type1 = factor(
            substr(colnames(counts)[5:100], 1, 2)
        ), Source = factor("ESC")
    )

    counts <- subset(counts,
        select = -c(
            EnsemblGeneID,
            EnsemblTranscriptID,
            AssociatedGeneName,
            GeneLength
        )
    )
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data, # look at split function
        rowData = row_data
    )
    counts <- assay(sce, "counts")
    libsizes <- colSums(counts)
    size.factors <- libsizes / mean(libsizes)
    logcounts(sce) <- log2(t(t(counts) / size.factors) + 1)
    return(sce)
}

g1_sce <- e_mtab("G1_singlecells_counts")
g2m_sce <- e_mtab("G2M_singlecells_counts")
s_sce <- e_mtab("S_singlecells_counts")
