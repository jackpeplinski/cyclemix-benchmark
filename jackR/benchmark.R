require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")

## QUESTIONS
# For e-mtab-3749, what is the format?
# For GSE, each file is the equivalent of a G1_cell#_count in E-MTAB?
# For marioni, is "ERCC-00004" ok to have in rows?

# don't need to log normalize, each file is a cell.
# Round fpkm for counts, and put fpkm into the logcounts
gse <- function() {
    files <- list.files(
        path =
            "./jackData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )
    for (file in files) {
        counts <- read.table(
            str_interp("./jackData/GSE42268/GSE42268_RAW/
            GSM1036480_EB5K_01.txt"),
            header = TRUE
        )

        rownames(counts) <- counts$id

        row_data <- DataFrame(
            row.names = counts$id,
            feature_symbol = factor(counts$gene.symbol)
        )

        counts <- subset(counts,
            select = -c(
                id,
                fpkm
            )
        )

        # col_data <- DataFrame(
        #     row.names = colnames(counts), # make my own G1_cell1_count?
        #     Species = factor("Mus musculus"),
        #     cell_type1 = factor("G1"), # all G1
        #     Source = factor("ESC") # still ESC?
        # )

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
}


emtab_2805 <- function(file_name) {
    counts <- read.table(
        str_interp("./jackData/E-MTAB-2805/E-MTAB-2805.processed.1/${file_name}.txt"),
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
        ),
        Source = factor("ESC")
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
    emtab_2805("G1_singlecells_counts"),
    emtab_2805("G2M_singlecells_counts"),
    emtab_2805("S_singlecells_counts")
)

# output2 <- classifyCells(e_mtab_res, MGeneSets$Cyclone)
# summary(factor(output2$phase))
# table(factor(output2$phase), e_mtab_res$cell_type1)
# plotMixture(output2$fit[["G2M"]], BIC = TRUE)

marioni <- function() {
    counts <- read.table(
        "./jackData/Marioni_lab_1_Jul_2015_gene_counts_table.txt"
    )

    rownames_length <- length(rownames(counts)) - 5
    row_data <- DataFrame(
        row.names = rownames(counts)[1:rownames_length]
        # feature_symbol = what to put here? should be a gene e.g., GNAI3
    )

    col_data <- DataFrame(
        row.names = colnames(counts),
        Species = factor("Mus musculus"), # check this
        cell_type1 = factor(), # what to put here? G1, S, G2M ...
        Source = factor("ESC")
    )

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data,
        rowData = row_data
    )
}
