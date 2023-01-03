require("SingleCellExperiment")
require("scater")
require("CycleMix")
require("stringr")

## Answered Questions
# For marioni, is "ERCC-00004" ok to have in rows? Don't need. RNA generated in the lab for control.

# Don't need to log normalize, each file is a cell.
# Round fpkm for counts, and put fpkm into the logcounts
gse <- function() {
    # counts should be something like
    # with headers: row.name (ensembleID),   gsmcellname
    #               e.g.,ENSMUSG00000000049  [fpkm]

    # load counts info; doesn't matter which file, could be any
    counts <- read.table(
        str_interp("./jackData/GSE42268/GSE42268_RAW/GSM1036480_EB5K_01.txt"),
        header = TRUE
    )

    # adds row.name with ensembleID
    rownames(counts) <- counts$id

    # adds gene.symbol with gene symbol
    row_data <- DataFrame(
        row.names = counts$id,
        feature_symbol = factor(counts$gene.symbol)
    )

    # remove unneeded info
    counts <- subset(counts,
        select = -c(
            id,
            gene.symbol,
            fpkm
        )
    )
    logcounts <- counts

    # copy/pasted this from gseParse
    col_data_xml <- data.frame(
        cell_type1 = c("G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "S", "S", "S", "S", "S", "S", "S", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1", "G1"), # nolint
        gsm = c("GSM1036483", "GSM1036484", "GSM1036485", "GSM1036486", "GSM1036487", "GSM1036488", "GSM1036489", "GSM1036490", "GSM1036491", "GSM1036492", "GSM1036493", "GSM1036494", "GSM1036495", "GSM1036496", "GSM1036500", "GSM1036501", "GSM1036502", "GSM1036503", "GSM1036504", "GSM1036505", "GSM1036506", "GSM1036507", "GSM1036508", "GSM1036509", "GSM1036510", "GSM1036511", "GSM1036512", "GSM1036513", "GSM1036522", "GSM1036523", "GSM1036524", "GSM1036525", "GSM1036526", "GSM1036527", "GSM1036528", "GSM1036529", "GSM1036530", "GSM1036531", "GSM1036532", "GSM1036533", "GSM1036534", "GSM1036535", "GSM1036536", "GSM1036537", "GSM1036538", "GSM1036539", "GSM1036540", "GSM1036541", "GSM1036542", "GSM1036543", "GSM1036544", "GSM1036545", "GSM1036546", "GSM1036547", "GSM1036548", "GSM1036549", "GSM1036550", "GSM1036551", "GSM1036552", "GSM1036553", "GSM1036554", "GSM1036555", "GSM1036556") # nolint
    )

    files <- list.files(
        path =
            "./jackData/GSE42268/GSE42268_RAW",
        pattern = ".*.txt"
    )

    for (file in files) {
        gsm <- substr(file, 1, 10)
        if (is.element(gsm, col_data_xml$gsm)) {
            counts_data <- read.table(
                str_interp("./jackData/GSE42268/GSE42268_RAW/${file}"),
                header = TRUE
            )

            format_counts <- function(count_type) {
                if (count_type != "logcounts") {
                    data <- data.frame(
                        lapply(
                            counts_data,
                            function(x) if (is.numeric(x)) round(x) else x
                        )
                    )
                } else {
                    data <- counts_data
                }

                rownames(data) <- data$id

                data <- subset(data,
                    select = -c(
                        id,
                        gene.symbol
                    )
                )

                colnames(data)[1] <- gsm
                return(data)
            }
            logcounts_data <- format_counts("logcounts")
            counts_data <- format_counts("counts")

            col_data <- DataFrame(
                row.names = col_data_xml$gsm,
                Species = factor("Mus musculus"),
                cell_type1 = col_data_xml$cell_type1, # why the "1" at the end?
                Source = factor("ESC") # this may be different
            )
            # print(col_data)
            # print(col_data["GSM1036531", ]) # verify that the correct phase is matched

            # print(head(counts_data))
            logcounts <- cbind(logcounts, logcounts_data)
            counts <- cbind(counts, counts_data)
        }
    }

    print(head(counts))
    # print(col_data)
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data, # look at split function
        rowData = row_data
    )
    counts <- assay(sce, "counts") # double check this and the logcounts is what she wanted.
    # print("ok")
    print(head(logcounts))
    logcounts(sce) <- logcounts
    return(sce)
}

gse_output <- gse()
output2 <- classifyCells(gse_output, MGeneSets$Cyclone)
summary(factor(output2$phase))
table(factor(output2$phase), gse_output$cell_type1)
plotMixture(output2$fit[["G2M"]], BIC = TRUE)

emtab_2805 <- function(file_name) {
    counts <- read.table(
        str_interp("./jackData/E-MTAB-2805/E-MTAB-2805.processed.1/${file_name}.txt"), # nolint
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
    # print(head(counts))
    counts <- subset(counts,
        select = -c(
            EnsemblGeneID,
            EnsemblTranscriptID,
            AssociatedGeneName,
            GeneLength
        )
    )
    print(head(counts))
    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data,
        rowData = row_data
    )
    counts <- assay(sce, "counts")
    libsizes <- colSums(counts)
    size.factors <- libsizes / mean(libsizes)
    logcounts(sce) <- log2(t(t(counts) / size.factors) + 1)
    # print(head(logcounts(sce)))
    return(sce)
}

emtab_2805_res <- cbind(
    emtab_2805("G1_singlecells_counts"),
    emtab_2805("G2M_singlecells_counts"),
    emtab_2805("S_singlecells_counts")
)

# output2 <- classifyCells(emtab_2805_res, MGeneSets$Cyclone)
# summary(factor(output2$phase))
# table(factor(output2$phase), emtab_2805_res$cell_type1)
# plotMixture(output2$fit[["G2M"]], BIC = TRUE)

marioni <- function() {
    counts <- read.table(
        "./jackData/Marioni_lab_1_Jul_2015_gene_counts_table.txt"
    )

    rownames_length <- length(rownames(counts)) - 5
    row_data <- DataFrame(
        row.names = rownames(counts)[1:rownames_length]
        # feature_symbol = what to put here. should be a gene e.g., GNAI3. Specify to use ENSMUSG00000000001 to  MGeneSets$Cyclone. Tried to use the example from EnsemblStuff. downloadEnsemblData() gives error, can't install biomaRt.
    )

    col_data <- DataFrame(
        row.names = colnames(counts),
        Species = factor("Mus musculus"),
        cell_type1 = factor(), # what to put here? G1, S, G2M ... look at the paper or database. Looked at both. It said just the three phases. There were two datasets mentioned...
        Source = factor("ESC")
    )

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = col_data,
        rowData = row_data
    )
}
