require("SingleCellExperiment")
require("scater")
require("CycleMix")

test <- read.table(
    "./jackData/E-MTAB-2805.processed.1/G1_singlecells_counts.txt",
    header = TRUE
)

counts <- read.table(
    "./jackData/E-MTAB-2805.processed.1/G1_singlecells_counts.txt",
    header = TRUE
)

rownames(counts) <- counts$EnsemblGeneID

counts <- head(counts, -97)
rowData <- DataFrame(
    row.names = counts$EnsemblGeneID,
    feature_symbol = factor(counts$AssociatedGeneName)
)

colData <- DataFrame(
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
    colData = colData, # look at split function
    rowData = rowData
)
counts <- assay(sce, "counts")
libsizes <- colSums(counts)
size.factors <- libsizes / mean(libsizes)
logcounts(sce) <- log2(t(t(counts) / size.factors) + 1)
